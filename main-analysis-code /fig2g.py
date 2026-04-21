import support_funcs as sf
import numpy as np
import matplotlib.pyplot as plt


def MULTIconditional_matching(data_dict,ses_n=1,plot=True):
    """
    Compute conditional matching metrics across mice for one session.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary of loaded data.
    ses_n : int, default=1
        Session number to analyze.
    plot : bool, default=True
        If True, plot summary conditional matching results.

    Returns
    -------
    None
    """
            
    # extracting just interval task, non opto mice 
    subset_dict = {
        subj: {'conc': data['conc']}
        for subj, data in data_dict.items()
        if 'conc' in data and not subj.startswith('6PO')
        }
    
    
    aport_match_rew,aport_match_vis = [],[]
    for i,mouse in enumerate(list(subset_dict.keys())):
        print(mouse)
        if data_dict[mouse]['conc']['b_meta']['nsessions'] >= ses_n:
            lick_df,_,_,bmeta,_                = sf.extract_data(data_dict,mouse,ses_n)  
            port_match_vis,port_match_rew      = conditional_matching(lick_df,bmeta,ses_n)
            aport_match_vis.append(port_match_vis), aport_match_rew.append(port_match_rew)
            
    if plot:
        plot_cond_match(aport_match_rew,aport_match_vis)


def conditional_matching(lick_df,bmeta,ses_n):
    """
    Compute conditional matching for each current port in one session.

    For each port, this function looks at the next port visited after leaving
    that port, then summarizes transition allocation in terms of both visit
    counts and rewarded counts. Matching is expressed as a ratio for each
    possible transition, conditioned on the current port.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Lick/visit dataframe for one session.
    bmeta : dict
        Behavioral metadata containing session interval information.
    ses_n : int
        Session number.

    Returns
    -------
    port_match_vis : numpy.ndarray
        Conditional matching ratios based on transition counts.
    port_match_rew : numpy.ndarray
        Conditional matching ratios based on rewarded transitions.

    Notes
    -----
    Transition matrices are sorted by the session's port intervals before
    being returned.
    """

    visits = lick_df[lick_df['unique_visit']==1][['port','rewarded','event_time']].reset_index(drop=True)

    port_match_vis,port_match_rew = [],[]
    for port in range(1,7):
        
        inext = np.array(visits[visits['port']==port].index+1)
        if inext[-1] == len(visits):
            inext = inext[:-1]
        next_vis = visits.iloc[inext]
        
        next_vis_count,next_rew_count = [],[]
        for i in np.sort(next_vis['port'].unique()): # for all the ports transitioned to from this port
            next_vis_count.append(len(next_vis[next_vis['port']==i])) # number of visits to that port (specific transition)
            next_rew_count.append(np.sum(next_vis[next_vis['port']==i]['rewarded'])) # number of rewards for that specific transition 
            
        # if not all transitions were made, insert NaNs - not zeros! 
        if len(list(set(list(range(1,7)))-set(np.sort(next_vis['port'].unique()))))>0: 
            imissing_trans = np.array(list((set(list(range(1,7)))-set(np.sort(next_vis['port'].unique())))))-1 # index of missing transitions
            for missingi in imissing_trans:
                next_vis_count = np.insert(np.array(next_vis_count).astype(float),missingi,np.nan)
                next_rew_count = np.insert(np.array(next_rew_count).astype(float),missingi,np.nan)
            
        # calculate matching in terms of ratio - but again only for specific transition away from the given port 
        pm_vis = np.array([next_vis_count[i]/(np.nansum(next_vis_count)-next_vis_count[i]) for i in range(len(next_vis_count))])
        pm_rew = np.array([next_rew_count[i]/(np.nansum(next_rew_count)-next_rew_count[i]) for i in range(len(next_rew_count))])
        
        pm_vis[pm_vis==float(0)] = 0.001
        pm_rew[pm_rew==float(0)] = 0.001
        
        port_match_vis.append(pm_vis)
        port_match_rew.append(pm_rew)
    
    # sort by interval
    port_match_vis = np.array(sf.sort_data_by_interval(port_match_vis,bmeta['intervals'][ses_n-1]))
    port_match_rew = np.array(sf.sort_data_by_interval(port_match_rew,bmeta['intervals'][ses_n-1]))
            
    return port_match_vis,port_match_rew


### PLOTTING

def plot_cond_match(aport_match_rew,aport_match_vis):
    """
    Plot conditional matching across mice on log-log axes.

    Each point represents a transition from one port to another, with
    reward-based matching on the x-axis and visit-based matching on the
    y-axis. Points are colored by the destination port.

    Parameters
    ----------
    aport_match_rew : list or array-like
        Per-mouse reward-based matching ratios (output of `conditional_matching`).
    aport_match_vis : list or array-like
        Per-mouse visit-based matching ratios.

    Returns
    -------
    None
    """
    color = plt.cm.plasma(np.linspace(0,0.9,6))
    hfont = {'fontname':'Helvetica'}

    fig, ax = plt.subplots(1,1,figsize=(3.8,3.5))

    for mouse in range(0,len(aport_match_rew)):
        for port in range(0,6):
            ax.scatter(np.log10(aport_match_rew[mouse][port]),np.log10(aport_match_vis[mouse][port]),color=color,s=10)
        ax.set_xlim([-2.5,1.5])
        ax.set_ylim([-2.5,1.5])
        ax.set_xticks([-2.5,-0.5,1.5])
        ax.set_yticks([-2.5,-0.5,1.5])
        ax.set_yticklabels([-2.5,-0.5,1.5])

    ax.plot(np.arange(-3.5,2),np.arange(-3.5,2),'gray',alpha=0.3,linewidth=2)
    subscript_10 = '\u2081\u2080'
    ax.set_xlabel('log' + subscript_10 + ' reward ratio',**hfont,fontsize='18')
    ax.set_ylabel('log' + subscript_10 + ' choice ratio',**hfont,fontsize='18')
    ax.spines[['right','top']].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.spines[['left','bottom']].set_linewidth(1.5)

    fig.tight_layout()
