import support_funcs as sf
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score,mean_squared_error


def eg_matching(data_dict,mouse='6PG6',ses_n=1):
    """
    Plot an example matching analysis for a single mouse and session.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary containing behavioral and licking data for all subjects.
    mouse : str, optional
        Subject identifier (default is '6PG6').
    ses_n : int, optional
        Session number to analyze (default is 1).

    Notes
    -----
    This function extracts licking and metadata for the specified mouse and
    session, then calls `m.matching_calc` with plotting enabled to visualize
    the matching relationship between reward and choice ratios.
    """
    
    lick_df  = data_dict[mouse]['conc']['mlick_df']
    bmeta    = data_dict[mouse]['conc']['b_meta']
    matching_calc(lick_df,bmeta,ses_n,plot_ind=True,data='real')



def MULTImatching(data_dict,ses_n=1,plot=True):
    """
    Compute matching behavior across multiple mice and optionally plot results.
    
    Parameters
    ----------
    data_dict : dict
        Nested dictionary containing behavioral datasets for all subjects.
    ses_n : int, optional
        Minimum session index required for inclusion (default is 1).
    plot : bool, optional
        If True, generate a summary matching plot across mice (default is True).
    
    Returns
    -------
    amatch_rew : list
        Reward matching values for each included mouse.
    amatch_vis : list
        Choice matching values for each included mouse.
    slopes : list
        Fitted slope values extracted from the matching analysis for each mouse.

    """
        
    # extracting just interval task, non opto mice 
    subset_dict = sf.subset_mice(data_dict, task='conc', include_opto=False, config=None, region=None)
    
    amatch_rew,amatch_vis,slopes = [],[],[]
    for mouse in subset_dict:
        print(mouse)
        if subset_dict[mouse]['conc']['b_meta']['nsessions'] >= ses_n:
            lick_df,_,_,bmeta,_                = sf.extract_data(data_dict,mouse,ses_n)  
            
            # standard matching and slopes
            match_rew,match_vis,fit_parameters = matching_calc(lick_df,bmeta,ses_n,plot_ind=False,data='real')
            if len(match_rew)>0:
                amatch_rew.append(match_rew), amatch_vis.append(match_vis)
                if len(fit_parameters)>0:
                    slopes.append(fit_parameters[0][0])
                else:
                    slopes.append(fit_parameters)
                    
    if plot:
        plot_match(amatch_rew,amatch_vis,data='real')
        
    return amatch_rew,amatch_vis,slopes 


def matching_calc(lick_df,bmeta,ses_n,plot_ind=False,data='real'):
    """
    Compute reward and choice matching ratios for a single session.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral dataframe containing visit information. For `data='real'`,
        this must include at least the columns `'unique_visit'`, `'port'`, and
        `'rewarded'`.
    bmeta : dict
        Session metadata dictionary containing interval information in
        `bmeta['intervals']`.
    ses_n : int
        Session number to analyze.
    plot_ind : bool, optional
        If True, plot the session-level matching results using
        `plot_ses_match` (default is False).
    data : {'real', 'AQUA'}, optional
        Type of input data. If `'real'`, visits are extracted from rows where
        `'unique_visit' == 1`. If `'AQUA'`, `lick_df` is assumed to already
        contain visit-level data (default is `'real'`).

    Returns
    -------
    match_rew : numpy.ndarray or list
        Reward ratio for each port, computed relative to the sum of rewards at
        all other ports. Returns an empty list if there are fewer than 100
        visits.
    match_vis : numpy.ndarray or list
        Visit ratio for each port, computed relative to the sum of visits at
        all other ports. Returns an empty list if there are fewer than 100
        visits.
    fit_parameters : list
        Matching fit summary in the form `[coef, r2, rmse]`, where `coef`
        contains the slope and intercept from the linear fit in log space.
        Returns an empty list if there are fewer than 100 visits.

    Notes
    -----
    Visits and rewards are first counted for each of the six ports, then sorted
    by session interval using `sf.sort_data_by_interval`. Matching is computed
    as the ratio for each port relative to all other ports combined. Zero and
    NaN values are replaced with 0.001 before log transformation. A linear fit
    is then performed on the log-transformed reward and visit ratios.
    """
    
    if data == 'real':
        visits      = lick_df[lick_df['unique_visit']==1][['port','rewarded']].reset_index(drop=True)
    elif data == 'AQUA':
        visits = lick_df 
    
    if len(visits)<100:
        match_rew,match_vis,fit_parameters = [],[],[]
        
    else:
        vis    = [len(visits[visits['port']==port]) for port in range(1,7)] # number of visits at that port 
        rew    = [visits[visits['port']==port]['rewarded'].cumsum().iloc[-1] for port in range(1,7)] # number of rewards at each port 
    
        # sort by interval
        vis = np.array(sf.sort_data_by_interval(vis,bmeta['intervals'][ses_n-1]))
        rew = np.array(sf.sort_data_by_interval(rew,bmeta['intervals'][ses_n-1]))

        match_vis = np.array([vis[i]/(np.sum(vis)-vis[i]) for i in range(len(vis))])
        match_rew = np.array([rew[i]/(np.sum(rew)-rew[i]) for i in range(len(rew))])
        
        # remove zeros or nans
        match_vis = np.nan_to_num(match_vis,nan=0.001)
        match_vis[match_vis==float(0)]=0.001 
        match_rew = np.nan_to_num(match_rew,nan=0.001)
        match_rew[match_rew==float(0)]=0.001 
        
        # fit line
        coef      = np.polyfit(np.log10(match_rew),np.log10(match_vis),1) # line of best fit: y = coef[0]*x + coef[1]. Least squares polynomial fit 
        poly1d_fn = np.poly1d(coef)
    
        # calculate r squared fit + RMSE
        ytrue = np.log10(match_vis) 
        ypred = coef[0]*np.log10(match_rew) + coef[1]
        r2    = r2_score(ytrue,ypred)
        rmse  = mean_squared_error(ytrue,ypred)
        
        fit_parameters = [coef,r2,rmse]
        
    if plot_ind:
        plot_ses_match(match_rew,match_vis,poly1d_fn,fit_parameters)
        
    return match_rew,match_vis,fit_parameters



### PLOTTING

def plot_match(amatch_rew,amatch_vis,data='real'):
    """
    Plot log10-transformed reward ratios vs. choice ratios across mice.
    
    Each mouse contributes multiple data points (e.g., transitions), plotted in log–log space.
    A unity line (y = x) is included to indicate perfect matching (choice ratio = reward ratio).
    
    Parameters
    ----------
    amatch_rew : list of array-like
        Reward ratios for each mouse (one array per subject). Values must be > 0.
    
    amatch_vis : list of array-like
        Choice (visit) ratios for each mouse (one array per subject). Values must be > 0.
    
    data : str, optional
        'real' or 'AQUA'. Controls axis limits and tick placement.
    
    Notes
    -----
    - Both axes are log10-transformed.
    - Points falling on the unity line reflect perfect matching.
    - Deviations from the line indicate under- or over-matching.
    """
    
    fig, ax = plt.subplots(1,1,figsize=(3.8,3.5))
    hfont = {'fontname':'Helvetica'}
    color = plt.cm.plasma(np.linspace(0,0.9,6))
    
    for mouse in range(0,len(amatch_rew)):
        ax.scatter(np.log10(amatch_rew[mouse]),np.log10(amatch_vis[mouse]),color=color,s=10)
        
    if data=='real':
        ax.set_xlim([-2,0.5])
        ax.set_ylim([-2,0.5])
        ax.set_xticks([-2.5,-1.5,-0.5,0.5])
        ax.set_yticks([-2.5,-1.5,-0.5,0.5])
        ax.set_yticklabels([-2.5,-1.5,-0.5,0.5])
    elif data=='AQUA':
        ax.set_xlim([-2,0.5])
        ax.set_ylim([-2,0.5])
        ax.set_xticks([-1.5,-0.5,0.5])
        ax.set_yticks([-1.5,-0.5,0.5])
        ax.set_yticklabels([-1.5,-0.5,0.5])
    
    ax.plot(np.arange(-2.0,1.0),np.arange(-2.0,1.0),'gray',alpha=0.3,linewidth=2)

    subscript_10 = '\u2081\u2080'
    ax.set_xlabel('log' + subscript_10 + ' reward ratio',**hfont,fontsize='18')
    ax.set_ylabel('log' + subscript_10 + ' choice ratio',**hfont,fontsize='18')
    ax.spines[['right','top']].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.spines[['left','bottom']].set_linewidth(1.5)

    fig.tight_layout()


def plot_ses_match(prew,pvis,poly1d_fn,fit_parameters):
    """
    Plot matching behavior for a single session with fitted and reference lines.

    Parameters
    ----------
    prew : array-like
        Reward ratios across conditions or transitions.
    pvis : array-like
        Choice (e.g., visual) ratios corresponding to `prew`.
    poly1d_fn : callable
        Fitted function (e.g., numpy.poly1d) used to generate the line of best fit
        in log10 space.
    fit_parameters : array-like
        Parameters of the fitted model (e.g., slope and intercept). Currently not
        displayed but included for completeness.

    Notes
    -----
    - Both reward and choice ratios are plotted in log10 space.
    - A dashed black line shows the fitted relationship.
    - A gray unity line (y = x) indicates perfect matching.
    """
    
    fig, ax = plt.subplots(1,1,figsize=(3.8,3.5))
    hfont = {'fontname':'Helvetica'}

    color = plt.cm.plasma(np.linspace(0,0.9,6))
    
    ax.plot(np.log10(prew),poly1d_fn(np.log10(prew)),'--k',linewidth=2,zorder=-1)
    ax.scatter(np.log10(prew),np.log10(pvis),color=color)
    x_range = np.log10([prew[0],prew[-1]])
    ax.plot(x_range,x_range,'gray',alpha=0.3,linewidth=2)
    ax.set_xlim([-2.5,0])
    ax.set_ylim([-2.5,0])
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




