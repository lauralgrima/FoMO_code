import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import support_funcs as sf
import scipy
import strategy_benchmarks as sb
import seaborn as sns
from scipy import stats

#----------------------------------------------------------------------------------
# illustrating deviation from random as a function of latency or n rewards
#----------------------------------------------------------------------------------

# functions to do with trying to calculate when (in visits or time) an animal shifts from random, including how random an animal is in the first place

def MULTIdiverge_time(data_dict,ses_n=1,nstrategy_samples=100,slope_thresh=0.02,win_len=10,ind_plot=True,plot=True):
    """
    Compute divergence-from-random metrics across mice for one session.

    This function subsets to non-opto mice with the 'conc' task, runs
    `diverge_rand()` for each mouse with the requested session available, and
    optionally plots divergence time and divergence reward number across mice.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary of loaded data.
    ses_n : int, default=1
        Session number to analyze.
    nstrategy_samples : int, default=100
        Number of random strategy simulations per mouse.
    slope_thresh : float, default=0.02
        Threshold on slope difference used to define divergence.
    win_len : int, default=10
        Moving window length, in visits.
    ind_plot : bool, default=True
        If True, plot per-mouse divergence results inside `diverge_rand()`.
    plot : bool, default=True
        If True, plot summary divergence metrics across mice.

    Returns
    -------
    diverge_lats_min : numpy.ndarray
        Divergence latencies, in minutes, for mice that diverged.
    i_no_diverge : list of int
        Indices of mice that did not diverge.
    """
    subset_dict = sf.subset_mice(data_dict, task='conc', include_opto=False, config=None, region=None)
    
    diverge_lats, diverge_lats_sem, diverge_viss, diverge_rews, diverge_rews_sem = [], [], [], [], []
    
    for mouse in subset_dict:
        if subset_dict[mouse]['conc']['b_meta']['nsessions'] >= ses_n:
            lick_df, _, _, bmeta, _ = sf.extract_data(data_dict, mouse, ses_n)
            diverge_lat, diverge_lat_sem, diverge_vis, diverge_rew, diverge_rew_sem, _ = diverge_rand(
                lick_df, ses_n, nstrategy_samples, bmeta,
                slope_thresh=slope_thresh, win_len=win_len,
                task='conc', plot=ind_plot
            )

            diverge_lats.append(diverge_lat)
            diverge_lats_sem.append(diverge_lat_sem)
            diverge_viss.append(diverge_vis)
            diverge_rews.append(diverge_rew)
            diverge_rews_sem.append(diverge_rew_sem)
    
    diverge_lats_min = np.array(list(filter(bool, diverge_lats))) / 60
    i_no_diverge     = [index for index, sublist in enumerate(diverge_lats) if not sublist]

    diverge_rews = np.array(list(filter(bool, diverge_rews)))

    if plot:
        plot_diverge(diverge_lats_min,'time from random (s)',[0,60,120,180],[-1,5],[0,200])
        plot_diverge(diverge_rews,'reward number from random',[0,50,100,150],[-1,5],[0,180])
            

    return diverge_lats_min,i_no_diverge


def diverge_rand(lick_df,ses_n,nstrategy_samples,bmeta,slope_thresh=0.02,win_len=10,task='conc',plot=False):
    """
    Estimate when behavior diverges from random reward collection.

    Divergence is defined as the first peak in the difference between the
    moving-window slopes of cumulative rewards for the animal and a simulated
    random strategy that exceeds `slope_thresh`.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Lick/visit dataframe for one session.
    ses_n : int
        Session number.
    nstrategy_samples : int
        Number of random strategy simulations to average over.
    bmeta : dict
        Behavioral metadata used for simulation.
    slope_thresh : float, default=0.02
        Threshold on slope difference used to define divergence.
    win_len : int, default=10
        Window length, in visits, for moving slope calculation.
    task : {'conc', 'prob'}, default='conc'
        Task type used to choose the random simulation.
    plot : bool, default=False
        If True, plot slope differences and cumulative rewards.

    Returns
    -------
    diverge_lat : float or list
        Latency to diverge from random, in seconds.
    diverge_lat_sem : float or list
        SEM of divergence latency.
    diverge_vis : list
        Visit indices at divergence.
    diverge_rew : float or list
        Number of cumulative rewards at divergence.
    diverge_rew_sem : float or list
        SEM of cumulative rewards at divergence.
    peak_df : pandas.DataFrame or list
        Peak slope differences and their times.

    Notes
    -----
    Returns empty lists if the animal does not outperform the simulated random
    strategy by the end of the session.
    """

    # mouse data 
    visit_times = np.array(lick_df.loc[lick_df['unique_visit']==1]['event_time'])
    cumu_rew    = lick_df.loc[lick_df['unique_visit']==1]['rewarded'].cumsum()
        
    # simulated random strategy
    rand_cumu_rew = []
    for n in range(nstrategy_samples+1):
        # random strategy
        print(n)
        if task == 'conc':
            rand_strategy_df = sb.simulate_strategies(visit_times,bmeta,ses_n,stagger=True,strategy='random')
        elif task == 'prob':
            rand_strategy_df = sb.simulate_strategies_prob(visit_times,bmeta,ses_n,strategy='random')
        rand_cumu_rew.append(rand_strategy_df['outcomes'].cumsum()) # cumulative random reward 

    rand_cumu_mean,rand_cumu_lCI,rand_cumu_hCI = sf.confidence_interval(rand_cumu_rew,confidence=0.95)
    
    # check that animal does end up with more reward than random. 
    if rand_cumu_mean[-1]>cumu_rew.iloc[-1]:
        diverge_lat,diverge_lat_sem,diverge_vis,diverge_rew,diverge_rew_sem,peak_df = [],[],[],[],[],[]

    else: 
        # fit slope to windows 
        win_anm   = sf.moving_window(win_len,cumu_rew)
        win_rand  = sf.moving_window(win_len,pd.Series(rand_cumu_mean))
        win_times = sf.moving_window(win_len,pd.Series(visit_times))
        
        anm_slopes  = np.array([np.polyfit(win_times[i], win_anm[i], 1)[0] for i in range(0,len(win_anm))])
        rand_slopes = np.array([np.polyfit(win_times[i], win_rand[i], 1)[0] for i in range(0,len(win_rand))])
    
        slope_diff = anm_slopes-rand_slopes
        
        # add initial value of 0 to help finding peaks
        slope_diff = np.insert(slope_diff,0,0)
        
        # find first peak on slope differences 
        peaks      = scipy.signal.find_peaks(slope_diff)[0]
        
        # create a df of slope differences at each peak, and include visit times 
        peak_slope_diffs = slope_diff[peaks]
        peak_times       = [np.mean(win_times[i]) for i in peaks]
        peak_df          = pd.DataFrame((peak_slope_diffs,peak_times)).transpose()
        peak_df.columns  = ['slope diffs','times']

        first_peak = peaks[np.where(slope_diff[peaks]>slope_thresh)[0]][0] # first peak over threshold
        
        # time to diverge from rand
        diverge_lat     = np.median(win_times[first_peak]) 
        diverge_lat_sem = stats.sem(win_times[first_peak])
        
        # number of visits to diverge from rand
        diverge_vis = [np.argmin(np.abs(visit_times-latency)) for latency in win_times[first_peak]]

        # how many rewards at divergence from rand 
        diverge_rew     = np.median(cumu_rew.iloc[diverge_vis])
        diverge_rew_sem = stats.sem(cumu_rew.iloc[diverge_vis])

        # individual plot (far left, fig 1e)
        if plot:
            if len(visit_times)>len(rand_cumu_mean):
                visit_times = visit_times[:-1]
                cumu_rew    = cumu_rew[:-1]
            
            plt.figure()
            plt.plot(slope_diff)
            plt.scatter(first_peak,slope_diff[first_peak])
            plt.scatter(peaks,slope_diff[peaks])
            
            fig,ax   = plt.subplots(figsize=(1.8,1.6))
            hfont  = {'fontname':'Arial'}
            ax.scatter(diverge_lat,diverge_rew,color='firebrick',zorder=1)
            ax.plot(visit_times,rand_cumu_mean,color='orange',zorder=-1)
            ax.plot(visit_times,cumu_rew,color='black',zorder=-1)
            ax.set_xticks([])
            ax.set_xticklabels([],fontsize=18,**hfont)
            ax.set_yticks([])
            ax.set_yticklabels([],fontsize=18,**hfont)
            ax.set_xlabel('time (min.)',**hfont,fontsize=18)
            ax.set_ylabel('cumu. rews',**hfont,fontsize=18)
            ax.spines[['right','top']].set_visible(False)
            ax.spines[['left','bottom']].set_linewidth(1.5)
            ax.set_xlim([0,10800])
            ax.set_ylim([0,500])
                
            fig.tight_layout()
            plt.show()

    return diverge_lat,diverge_lat_sem,diverge_vis,diverge_rew,diverge_rew_sem,peak_df


def plot_diverge(diverge_metric,ylabel,yticks,xlim,ylim):
    """
    Plot the distribution of a single divergence metric as a boxplot with
    overlaid individual data points.

    This function displays a summary of divergence-from-random values using a
    boxplot (without fliers or caps) together with a stripplot showing the
    individual observations. 

    Parameters
    ----------
    diverge_metric : array-like
        Values to plot.
    ylabel : str
        Y-axis label.
    yticks : array-like
        Y-axis tick locations.
    xlim : tuple
        X-axis limits.
    ylim : tuple
        Y-axis limits.


    Returns
    -------
    None
    
    """
    
    hfont    = {'fontname':'Arial'} 
    fig,ax   = plt.subplots(figsize=(3,3))

    dot_colors  = len(diverge_metric)*['silver']

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    sns.boxplot(y=diverge_metric,x=0,ax=ax,color='black',fill=False,width=0.5,showfliers=False,showcaps=False)
    sns.stripplot(y=diverge_metric,x=0,edgecolor=dot_colors,facecolor='white',linewidth=2,zorder=-1)
    
    ax.set_ylabel(ylabel,**hfont,fontsize=18)
    
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=18)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)

    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()
    plt.show()


        
        


