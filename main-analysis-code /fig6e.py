import support_funcs as sf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pypalettes import load_cmap


def MULTIses_rewpeak_DA(data_dict,tw_start=0,tw_length=2,n_quants=4):
    '''
    Compute reward-evoked dopamine responses across sessions and animals.

    For each mouse and session, extracts photometry responses to reward events,
    computes mean response magnitude and quantile-based summaries, and tracks
    reward-evoked signals over time. Signals are optionally smoothed and
    interpolated to create continuous time series for each animal. Results are
    aggregated across mice (NAc and DMS separated) and used for visualization
    and downstream analyses (e.g., relating DA responses to behavior or learning).

    Parameters
    ----------
    data_dict : dict
        Nested data structure containing behavioral, photometry, and metadata
        for all mice and sessions.
    tw_start : float, optional
        Start of the time window (in seconds) relative to reward used for
        extracting responses. Default is 0.
    tw_length : float, optional
        Length of the time window (in seconds) for response extraction.
        Default is 2.
    n_quants : int, optional
        Number of quantiles used to group reward responses within each session.
        Default is 4.

    Returns
    -------
    NAc_ind_rews : list of pandas.DataFrame
        Time-aligned, interpolated reward-evoked dopamine signals for each NAc mouse.
    NAc_quart_sig : list of numpy.ndarray
        Quantile-binned reward response magnitudes for each NAc mouse across sessions.

    Notes
    -----
    - Mean and quantile response metrics are computed per session and then aggregated.
    - Reward responses are smoothed (rolling window) and interpolated to produce
      continuous time series.
    - Only mice with more than two sessions are included.
    - Plotting is performed internally via `plot_multises_rewpeak_DA`.
    '''
    subset_dict = sf.subset_mice(data_dict,config=1, region='NAc')

    NAc_mean_sig, NAc_quart_sig, NAc_ind_rews = [],[],[]
    DMS_mean_sig, DMS_quart_sig, DMS_ind_rews = [],[],[]
    for mouse in list(subset_dict.keys()):
        anm_mean_sig,anm_quart_sig,anm_ind_rews = [],[],[]
        print(mouse)
        
        tot_ses   = data_dict[mouse]['conc']['b_meta']['nsessions'] # total number of sessions for given mouse 
        ses_count = list(range(1,tot_ses+1))
        if len(ses_count)>2: # if animal did more than just the initial two sessions
            for ses in ses_count:
                lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,ses)
                mean_sig, quart_sig,ind_rews        = ses_rewpeak_DA(lick_df,photo_df,pmeta,bmeta,ses,tw_start,tw_length,n_quants,peak_or_mean='mean',plot=False)
                ind_rews_smooth                     = pd.DataFrame(ind_rews['resp'].rolling(window=31).mean()) # smooth for specific session
                ind_rews_smooth['event_time']       = ind_rews['event_time']
                bin_edges  = np.arange(0,10800+1,1)
                bin_counts = np.histogram(ind_rews_smooth['event_time'],bins=bin_edges)[0]
                # interpolate by time 
                resp_count = 0
                long_events = []
                for i,count in enumerate(bin_counts):
                    if bin_counts[i]==1:
                        long_events.append(ind_rews_smooth['resp'].iloc[resp_count])
                        resp_count+=1
                    else:
                        long_events.append(np.nan)
                ind_rews_interpolate = pd.DataFrame(long_events).interpolate(method='linear')
                anm_mean_sig.append(mean_sig),anm_quart_sig.append(quart_sig),anm_ind_rews.append(ind_rews_interpolate)
            # join dfs across sessions and interpolate
            anm_ind_rews = pd.concat(anm_ind_rews).interpolate(method='linear').reset_index(drop=True)
            if pmeta['region'] == 'NAc':
                NAc_mean_sig.append(np.array(anm_mean_sig))
                NAc_quart_sig.append(np.array(anm_quart_sig))
                DMS_mean_sig.append(np.array(anm_mean_sig))
                NAc_ind_rews.append(anm_ind_rews)
            elif pmeta['region'] == 'DMS':
                DMS_quart_sig.append(np.array(anm_quart_sig))
                DMS_ind_rews.append(anm_ind_rews)
        elif len(ses_count)<2:
            continue 

    # # reorganise
    NAc_ms = np.array(pd.DataFrame(NAc_mean_sig).transpose())
    NAc_qs = [np.concatenate(NAc_quart_sig[i]) for i in range(len(NAc_quart_sig))]
    NAc_qs = np.array(pd.DataFrame(NAc_qs).transpose())
    
    # calculate median for plotting individual animals 
    NAc_ind_median = pd.concat(NAc_ind_rews,axis=1).median(axis=1)

    # stats on quartiles 
    longform_qs         = pd.melt(pd.DataFrame(NAc_qs).transpose())
    longform_qs['subj'] = np.concatenate([np.arange(1,len(NAc_qs[0])+1)]*20)
    longform_qs['ses']  = np.concatenate([np.repeat(i,len(NAc_qs[0])*4) for i in [1,2,3,4,5]])
    
    # plotting
    plot_multises_rewpeak_DA(NAc_ms,NAc_qs,NAc_ind_rews,NAc_ind_median)
    
    return NAc_ind_rews, NAc_quart_sig

def ses_rewpeak_DA(lick_df,photo_df,pmeta,bmeta,n_ses,tw_start=0,tw_length=2,n_quants=4,peak_or_mean='mean',plot=False):
    '''
    Compute reward-evoked dopamine response magnitudes for a single session.

    Extracts photometry signal windows aligned to rewarded unique visits and
    summarizes each response as either the peak or mean signal within the specified
    time window. Responses are also averaged within sequential quantiles of reward
    events to assess changes across the session.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral event data containing reward outcome, unique visit labels,
        event times, port identity, and photometry sample indices.
    photo_df : pandas.DataFrame
        Photometry data containing the z-scored dopamine signal.
    pmeta : dict
        Photometry metadata containing sampling rate information for each session.
    bmeta : dict
        Behavioral metadata, used for optional plot titles.
    n_ses : int
        Session number used to index the corresponding sampling rate.
    tw_start : float, optional
        Start of the response window, in seconds relative to reward/visit.
        Default is 0.
    tw_length : float, optional
        Length of the response window, in seconds. Default is 2.
    n_quants : int, optional
        Number of sequential quantiles used to summarize reward responses
        across the session. Default is 4.
    peak_or_mean : {'peak','mean'}, optional
        Whether to summarize each reward response by its peak or mean signal
        within the response window. Default is 'mean'.
    plot : bool, optional
        If True, plots individual reward response magnitudes for the session.
        Default is False.

    Returns
    -------
    mean_sig : float
        Mean reward-evoked response magnitude across all rewarded visits.
    quart_sig : list of float
        Mean reward-evoked response magnitude within each sequential quantile
        of rewarded visits.
    ind_rews : pandas.DataFrame
        Reward-level response data containing event time and response magnitude
        for each rewarded visit.
    '''
    photo          = photo_df['signal'].to_numpy()
    photo_i_rew    = lick_df.loc[(lick_df['rewarded']==1) & (lick_df['unique_visit']==1)][['photo_i','port']].dropna()
    win_idx_rew    = sf.extract_window(photo_i_rew['photo_i'],photo,pmeta['sampling_rate'][n_ses-1],tw_start=tw_start,tw_length=tw_length)[0]
    win_signal_rew = np.vstack([photo[win_idx_rew[x][0]:win_idx_rew[x][1]] for x in range(len(win_idx_rew))])
    
    if peak_or_mean == 'peak':
        max_signal_rew = [np.max(win_signal_rew[i]) for i in range(len(win_signal_rew))]
    elif peak_or_mean == 'mean':
        max_signal_rew = [np.mean(win_signal_rew[i]) for i in range(len(win_signal_rew))]
    mean_sig = np.mean(max_signal_rew)
    
    if plot: 
        plt.figure()
        plt.title(bmeta['animal_no']+' ' +'ses' + ' ' + str(n_ses))
        plt.scatter(range(len(max_signal_rew)),max_signal_rew)
        
    # splitting into quantiles 
    bin_edges      = np.linspace(0, len(max_signal_rew), n_quants+1).astype(int)
    quart_sig      = [np.mean(max_signal_rew[bin_edges[i]:bin_edges[i+1]]) for i in range(len(bin_edges)-1)]
    
    # returning each individual reward response, plus event time for plotting 
    ind_rews         = pd.DataFrame(lick_df[(lick_df['rewarded']==1)&(lick_df['unique_visit']==1)][['event_time','photo_i']]).dropna()
    ind_rews         = pd.DataFrame(ind_rews['event_time'])
    ind_rews['resp'] = max_signal_rew
    
    return mean_sig, quart_sig, ind_rews
    


def plot_multises_rewpeak_DA(mean_sig,quart_sig,ind_anm,ind_median):
    '''
    Plot reward-evoked dopamine responses across sessions.

    Generates three summary plots: session-level mean reward responses,
    within-session quantile-binned reward responses, and continuous reward-response
    trajectories across sessions. Individual mice are shown in gray and group-level
    summaries are shown using boxplots or a median trace.

    Parameters
    ----------
    mean_sig : array-like
        Session-level mean reward response magnitudes. Each row/session contains
        one value per mouse.
    quart_sig : array-like
        Quantile-binned reward response magnitudes across sessions. Each entry
        contains response values for one session quantile across mice.
    ind_anm : list of pandas.DataFrame
        Interpolated reward-response trajectories for individual mice.
    ind_median : pandas.Series
        Median interpolated reward-response trajectory across mice.

    Returns
    -------
    None
        Generates matplotlib figures.
    '''
    hfont = {'fontname':'Arial'}
    fig,ax = plt.subplots(figsize=(3,3))

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    cmap   = load_cmap("Starfish")
    mcolors = [cmap(i) for i in [0,0.25,0.5,0.75,1.0]]
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    xpairs = [[0,0.1],[0.2,0.3],[0.4,0.5],[0.6,0.7],[0.8,0.9]]
    for i in range(len(mean_sig)):
        sns.boxplot(y=mean_sig[i],x=xpairs[i][0],ax=ax,color=mcolors[i],fill=False,width=0.4,showfliers=False,showcaps=False)
     #   sns.stripplot(y=mean_sig[i],x=xpairs[i][0],edgecolor='gainsboro',facecolor='white',linewidth=2,zorder=-1)
        
    [sns.lineplot(x=[0,1,2,3,4],y=[mean_sig[0][i],mean_sig[1][i],mean_sig[2][i],mean_sig[3][i],mean_sig[4][i]],ax=ax,color='gainsboro',linewidth=1.5,zorder=-1) for i in range(0,len(mean_sig[0]))]
    
    ax.set_ylabel('mean DA to reward \n(z-score)',**hfont,fontsize=18)
    
    ax.set_yticks([0,0.4,0.8,1.2])
    ax.set_yticklabels([0,0.4,0.8,1.2],fontsize=18)
    ax.set_xlim([-1,5])
    ax.set_ylim([0,1.2])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()

    fig,ax = plt.subplots(figsize=(8,3))
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    xpairs  = np.arange(0,4,0.2)
    qcolors = [cmap(i) for i in np.concatenate([[i]*4 for i in [0,0.25,0.5,0.75,1]])]
    
    for j in range(len(quart_sig[0])):
        ys = [quart_sig[i][j] for i in range(len(quart_sig))]
        [sns.lineplot(x=list(range(0,20)),y=ys,ax=ax,color='gainsboro',linewidth=1.5)]

    for i in range(len(quart_sig)):
        sns.boxplot(y=quart_sig[i],x=xpairs[i],ax=ax,color=qcolors[i],fill=False,width=0.4,showfliers=False,showcaps=False)
       # sns.stripplot(y=quart_sig[i],x=xpairs[i],edgecolor='gainsboro',facecolor='white',linewidth=1.5,zorder=-1)

    ax.set_ylabel('mean DA to reward \n(z-score)',**hfont,fontsize=18)
    
    ax.set_yticks([0,0.4,0.8,1.2])
    ax.set_yticklabels([0,0.4,0.8,1.2],fontsize=18)
    ax.set_xlim([-1,20])
    ax.set_ylim([-0.2,1.4])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()
    
    fig,ax1 = plt.subplots(figsize=(7,3))
    ax2     = ax1.twiny()
    
    # adding quartiles
    edge_wins = np.arange(0,54000,2700)
    q_x       = (edge_wins[:-1]+edge_wins[1:])/2
    q_x       = np.append(q_x,q_x[1]-q_x[0]+q_x[-1])
    
    [sns.lineplot(x=np.arange(0,len(ind_anm[i]),1),y=np.array(ind_anm[i][0].rolling(window=100,min_periods=10).mean()),linewidth=0.8,color='gainsboro',ax=ax1,zorder=-1) for i in range(len(ind_anm))]
    sns.lineplot(x=np.arange(0,len(ind_median)),y=ind_median.rolling(window=500,min_periods=10).mean(),color='black',linewidth=2,ax=ax1)
    
  #  fig,ax2 = plt.subplots(figsize=(7,3))
    for i in range(len(quart_sig)):
        sns.boxplot(y=quart_sig[i],x=i+1,color=qcolors[i],fill=False,width=0.4,showfliers=False,showcaps=False,ax=ax2,zorder=10,linewidth=1.5) 
        
    ax2.set_xlim([-1.8,21])
    ax1.set_ylim([-0.3,1.8])
    ax1.set_yticks([-0.3,1.8])
    ax1.set_yticklabels([-0.3,1.8],fontsize=18)
    ax1.set_xticklabels([])
    ax2.spines[['top','right','bottom']].set_visible(False)
    ax1.spines[['top','right','bottom']].set_visible(False)
    ax2.set_xticks([])
    ax1.set_xticks([])
    ax1.spines['left'].set_linewidth(1.5)
    
    cmap    = load_cmap("Starfish")
    sescolors = [cmap(i) for i in [0.25,0.5,0.75,1.0]]
    ax2.axvline(x=3.5, ymin=0,ymax=1,linestyle='dotted',color=sescolors[0])
    ax2.axvline(x=7.5,ymin=0,ymax=1,linestyle='dotted',color=sescolors[1])
    ax2.axvline(x=11.5,ymin=0,ymax=1,linestyle='dotted',color=sescolors[2])
    ax2.axvline(x=15.5,ymin=0,ymax=1,linestyle='dotted',color=sescolors[3])

    fig.tight_layout()
    

