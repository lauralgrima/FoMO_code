import numpy as np
import support_funcs as sf
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score,mean_squared_error


def eg_learning_to_match(data_dict,mouse='6PG6',ses_n=1):
    """
    Convenience wrapper to plot learning-to-match for a single mouse.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary containing behavioral data for all subjects.
    mouse : str
        Subject identifier.
    ses_n : int
        Session number to analyze.

    """
    
    lick_df  = data_dict[mouse]['conc']['mlick_df']
    bmeta    = data_dict[mouse]['conc']['b_meta']

    learning_to_match(lick_df,bmeta,ses_n,win_len=45,plot_ind=True)
    
def MULTIlearning_to_match(data_dict,ses_n=1,win_len=30,plot=True):
    """
    Compute windowed matching sensitivity across multiple mice.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary containing behavioral data for all subjects.
    ses_n : int
        Session number to analyze.
    win_len : int, optional
        Length of each time window in minutes used for learning-to-match
        analysis (default is 30).
    plot : bool, optional
        If True, plot sensitivity across windows for all included mice
        (default is True).

    Returns
    -------
    None
    """
    
    # extracting just interval task, non opto mice 
    subset_dict = sf.subset_mice(data_dict, task='conc', include_opto=False, config=None, region=None)

    win_sensitivity = []
    for i,mouse in enumerate(list(subset_dict.keys())):
        print(mouse)
        if subset_dict[mouse]['conc']['b_meta']['nsessions'] >= ses_n:
            lick_df,_,_,bmeta,_ = sf.extract_data(data_dict,mouse,ses_n)  

            # # windowed matching and slopes 
            _,_,win_fit_params = learning_to_match(lick_df,bmeta,ses_n,win_len=win_len,plot_ind=False)
            win_sensitivity.append(np.array([window[0][0] for window in win_fit_params]))
            
    if plot:
        plot_win_sensitivity(np.vstack(win_sensitivity))
        

def learning_to_match(lick_df, bmeta, ses_n, win_len=30, plot_ind=False):
    """
    Compute the evolution of matching sensitivity across a session in fixed time windows.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral dataframe containing visit information. Must include
        `'unique_visit'`, `'port'`, `'rewarded'`, and `'event_time'` columns.
    bmeta : dict
        Session metadata dictionary containing interval information in
        `bmeta['intervals']`.
    ses_n : int
        Session number to analyze.
    win_len : int, optional
        Length of each time window in minutes (default is 30).
    plot_ind : bool, optional
        If True, plot matching for each time window using
        `plot_ses_match_learn` (default is False).

    Returns
    -------
    win_match_vis : list of numpy.ndarray
        Visit matching ratios for each time window.
    win_match_rew : list of numpy.ndarray
        Reward matching ratios for each time window.
    win_fit : list of list
        Fit parameters for each time window in the form `[coef, r2, rmse]`,
        where `coef` contains the slope and intercept of the fitted line in
        log10 space.

    Notes
    -----
    The session is divided into consecutive fixed windows of length `win_len`.
    Within each window, visits and rewards are counted across the six ports,
    sorted by session interval, and converted into matching ratios relative to
    the sum across the other ports. A linear fit is then computed in log10
    space, and the slope of this fit can be used as a measure of matching
    sensitivity over time.
    """
    visits = lick_df[lick_df['unique_visit'] == 1][['port', 'rewarded', 'event_time']].reset_index(drop=True)

    windows = [i * win_len * 60 for i in range(0, int(10800 / (win_len * 60)))]
    windows.append(windows[-1] + (win_len * 60))

    win_visits = [
        visits[(visits['event_time'] > windows[i]) & (visits['event_time'] <= windows[i + 1])]
        for i in range(0, len(windows) - 1)
    ]

    win_match_vis, win_match_rew, win_fit = [], [], []
    poly1d_fns, afit_parameters = [], []

    for win in win_visits:
        vis = [len(win[win['port'] == port]) for port in range(1, 7)]
        rew = []

        for port in range(1, 7):
            if not win[win['port'] == port]['rewarded'].cumsum().empty:
                port_rew = win[win['port'] == port]['rewarded'].cumsum().iloc[-1]
            else:
                port_rew = 0.0
            rew.append(port_rew)

        vis = np.array(sf.sort_data_by_interval(vis, bmeta['intervals'][ses_n - 1]), dtype=float)
        rew = np.array(sf.sort_data_by_interval(rew, bmeta['intervals'][ses_n - 1]), dtype=float)
        
        match_vis = np.array([vis[i]/(np.sum(vis)-vis[i]) for i in range(len(vis))])
        match_rew = np.array([rew[i]/(np.sum(rew)-rew[i]) for i in range(len(rew))])
        
        # remove zeros or nans
        match_vis = np.nan_to_num(match_vis,nan=0.001)
        match_vis[match_vis==float(0)]=0.001 
        match_rew = np.nan_to_num(match_rew,nan=0.001) 
        match_rew[match_rew==float(0)]=0.001 

        win_match_vis.append(match_vis)
        win_match_rew.append(match_rew)

        x = np.log10(match_rew)
        y = np.log10(match_vis)

        coef = np.polyfit(x, y, 1)
        poly1d_fn = np.poly1d(coef)

        ypred = poly1d_fn(x)
        r2 = r2_score(y, ypred)
        rmse = mean_squared_error(y, ypred)

        fit_parameters = [coef, r2, rmse]

        win_fit.append(fit_parameters)
        poly1d_fns.append(poly1d_fn)
        afit_parameters.append(fit_parameters)

    if plot_ind:
        plot_ses_match_learn(win_match_rew, win_match_vis, poly1d_fns, afit_parameters)

    return win_match_vis, win_match_rew, win_fit



### PLOTTING

def plot_ses_match_learn(win_match_rew,win_match_vis,poly1d_fns,afit_parameters):
    """
    Plot fitted matching relationships across time windows within a session.

    Parameters
    ----------
    win_match_rew : list of numpy.ndarray
        Reward matching ratios for each time window.
    win_match_vis : list of numpy.ndarray
        Visit matching ratios for each time window.
    poly1d_fns : list of callable
        Fitted linear functions for each time window, typically generated with
        `numpy.poly1d` in log10 space.
    afit_parameters : list of list
        Fit summaries for each time window in the form `[coef, r2, rmse]`.
        Currently included for completeness but not directly displayed.

    Notes
    -----
    This function overlays the fitted matching lines for multiple time windows
    from a single session on the same axes. Each window is plotted in a
    different shade of green, and a gray unity line (y = x) is included as a
    reference for perfect matching. The axes are shown in log10 space.
    """
    
    fig, ax = plt.subplots(1,1,figsize=(3.8,3.5))
    hfont = {'fontname':'Helvetica'}
    line_color = plt.cm.Greens(np.linspace(0.3,0.8,len(win_match_rew)))
    
    for win in range(0,len(win_match_rew)):
       ax.plot(np.log10(win_match_rew[win]),poly1d_fns[win](np.log10(win_match_rew[win])),'--',color=line_color[win],linewidth=2)
       
    x_range = np.log10([win_match_rew[0][0],win_match_rew[0][-1]])
    ax.plot(x_range,x_range,'gray',alpha=0.3)
       
    ax.set_xlim([-2.5,0.5])
    ax.set_ylim([-2.5,0.5])
    
    ax.set_xticks([-2.5,-1.5,-0.5,0.5])
    ax.set_yticks([-2.5,-1.5,-0.5,0.5])
    ax.set_yticklabels([-2.5,-1.5,-0.5,0.5])
    subscript_10 = '\u2081\u2080'
    ax.set_xlabel('log' + subscript_10 + ' reward ratio',**hfont,fontsize='18')
    ax.set_ylabel('log' + subscript_10 + ' choice ratio',**hfont,fontsize='18')
    ax.spines[['right','top']].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.spines[['left','bottom']].set_linewidth(1.5)

    fig.tight_layout()


def plot_win_sensitivity(win_sensitivity):
    """
    Plot the evolution of matching sensitivity across time windows.

    Parameters
    ----------
    win_sensitivity : array-like
        2D array or list of shape (n_subjects, n_windows) containing
        sensitivity (e.g., slope) values for each time window and subject.

    Notes
    -----
    - Displays a boxplot for each time window showing the distribution
      of sensitivity values across subjects.
    - Overlays the mean sensitivity across subjects as a line with points.
    - Time windows are assumed to be evenly spaced and are labeled in minutes
      (e.g., 30, 60, ..., 180).
    """
    
    mean_sensitivity = np.mean(win_sensitivity,axis=0)
        
    fig, ax    = plt.subplots(1,1,figsize=(6,3.5))
    hfont      = {'fontname':'Helvetica'}

    win_sens   = [[win_sensitivity[i][j] for i in range(0,len(win_sensitivity))] for j in range(0,len(win_sensitivity[0]))]
    line_color = plt.cm.Greens(np.linspace(0.3,0.8,len(win_sens)))
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    ax.plot(mean_sensitivity,color='dimgrey',zorder=-1)
    ax.scatter([0,1,2,3,4,5],mean_sensitivity,color=line_color,s=50)
    for quartile in range(0,len(win_sens)):
        sns.boxplot(y=win_sens[quartile],x=quartile,ax=ax,native_scale=True,width=0.2,fill=False,color=line_color[quartile],showfliers=False,showcaps=False)

    ax.set_ylim([0,1])
    ax.set_xlim([-0.5,len(win_sens)])
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_yticklabels([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_xticks([0,1,2,3,4,5])
    ax.set_xticklabels([30,60,90,120,150,180])
    ax.set_ylabel('sensitivity',**hfont,fontsize=18)
    ax.set_xlabel('time in session (min.)',**hfont,fontsize=18)

    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')

    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))

    fig.tight_layout()

