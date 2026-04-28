import support_funcs as sf
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score,root_mean_squared_error


def MULTImatch_opto_4worst(data_dict,limit=90,ses_n=1,plot=True):
    """
    Compare matching behavior between optogenetic and control mice.
    
    Computes matching metrics (reward-based and visit-based) for opto and control
    groups, optionally restricting analysis to behavior after a time threshold.
    Designed for “worst 4” opto animals (excluding specified non-naive mice).
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral dataset.
    limit : float or None, default=90
        Time threshold in minutes. Only visits occurring after this time are
        included. If None or 0, all visits are used.
    ses_n : int, default=1
        Session number to analyze.
    plot : bool, default=True
        Whether to plot matching results across groups.
    
    Returns
    -------
    None
    """
    
    opto_dict  = sf.subset_mice(data_dict,opto_only=True)
    opto_dict  = {k: v for k, v in opto_dict.items() if k not in ['6PO3','6PO4']} # remove non-naive mice
    cntrl_dict = sf.subset_mice(data_dict,config=1,task='conc')

    ematch_rew,ematch_vis = [],[]
    for exp in [opto_dict,cntrl_dict]:
        mmatch_rew,mmatch_vis = [],[]
        for mouse in exp:
            print(mouse)
            lick_df,_,_,bmeta,_ = sf.extract_data(data_dict,mouse,ses_n)  # session 1 only
            visits              = lick_df[lick_df['unique_visit']==1][['port','rewarded','event_time','session_number']].reset_index(drop=True)
            if limit:
                visits = visits[visits['event_time']>limit*60][['port','rewarded']].reset_index(drop=True)
            match_rew,match_vis,fit_parameters = simple_matching(visits,bmeta,ses_n)
            mmatch_rew.append(match_rew)
            mmatch_vis.append(match_vis)
        ematch_rew.append(mmatch_rew)
        ematch_vis.append(mmatch_vis)
        
    if plot:
        plot_all_match(ematch_rew,ematch_vis)


def simple_matching(visits,bmeta,ses_n):
    """
    Compute simple matching between reward and visit ratios across ports.
    
    Calculates, for each port, the ratio of rewards earned and visits made at
    that port relative to all other ports. Optionally sorts ports by interval,
    then fits a log-log matching relationship between reward ratio and visit ratio.
    
    Parameters
    ----------
    visits : pandas.DataFrame
        Visit-level data containing 'port' and 'rewarded' columns.
    bmeta : dict-like or None
        Behavioral metadata containing session intervals. If provided, ports are
        sorted by interval.
    ses_n : int
        Session number.
    
    Returns
    -------
    tuple
        Reward ratios, visit ratios, and fit parameters:
        [linear fit coefficients, R², RMSE].
    """
    vis    = [len(visits[visits['port']==port]) for port in range(1,7)] # number of visits at that port 
    rew    = visits.groupby("port")["rewarded"].sum()
    rew    = rew.reindex(range(1, 7), fill_value=0).tolist()

    # sort by interval
    if bmeta:
        if bmeta['intervals'][ses_n-1][0]%1==0:
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

    # calculate r squared fit + RMSE
    ytrue = np.log(match_vis) 
    ypred = coef[0]*np.log(match_rew) + coef[1]
    r2    = r2_score(ytrue,ypred)
    rmse  = root_mean_squared_error(ytrue,ypred)
    
    fit_parameters = [coef,r2,rmse]

    return match_rew,match_vis,fit_parameters

# PLOTTING

def plot_all_match(ematch_rew,ematch_vis):
    """
    Plot matching behavior across animals for experimental groups.
    
    Displays log-transformed reward ratios versus visit ratios for each port
    (6 points per animal), comparing two groups (e.g., opto vs control). A unity
    line is included for reference.
    
    Parameters
    ----------
    ematch_rew : list
        Nested list of reward ratios. Outer list indexes experimental groups,
        inner list indexes animals, and each element contains an array of length 6.
    ematch_vis : list
        Nested list of visit ratios with the same structure as `ematch_rew`.
    
    Returns
    -------
    None
    """
    
    fig, ax = plt.subplots(1,1,figsize=(3.8,3.5))
    hfont = {'fontname':'Helvetica'}
    
    # for real data, -2.5 to 1.0
    ax.plot(np.arange(-2.5,1.0),np.arange(-2.5,1.0),'gray',alpha=0.3,linewidth=2,zorder=-10)

    for exp in range(len(ematch_rew)):
        amatch_rew = ematch_rew[exp]
        amatch_vis = ematch_vis[exp]
        
        for mouse in range(0,len(amatch_rew)):
            if exp == 0:
                ax.scatter(np.log10(amatch_rew[mouse]),np.log10(amatch_vis[mouse]),color=['navy','navy','deepskyblue','deepskyblue','deepskyblue','deepskyblue'],s=15,zorder=10)    
            elif exp == 1:
                ax.scatter(np.log10(amatch_rew[mouse]),np.log10(amatch_vis[mouse]),color='lightgrey',s=15,zorder=-10)
    ax.set_xlim([-2.5,0.5])
    ax.set_ylim([-2.5,0.5])
    ax.set_xticks([-2.5,0.5])
    ax.set_yticks([-2.5,0.5])
    ax.set_yticklabels([-2.5,0.5])

    ax.set_xlabel('log$_{10}$ reward ratio',**hfont,fontsize='18')
    ax.set_ylabel('log$_{10}$ choice ratio',**hfont,fontsize='18')
    ax.spines[['right','top']].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.spines[['left','bottom']].set_linewidth(1.5)
        
    fig.tight_layout()
    plt.show()
