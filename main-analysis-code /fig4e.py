import support_funcs as sf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def eg_ind_resp(data_dict,mouse='6PG6',tw_start=0,tw_length=2,ses_n=1,plot=True):
    """
    Plot example event-level dopamine responses for one mouse.
    
    Extracts rewarded and unrewarded port-specific dopamine response peaks for a
    single session and optionally plots individual responses with smoothed traces.
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral and photometry data.
    mouse : str, default='6PG6'
        Mouse ID to plot.
    tw_start : float, default=0
        Window start time relative to visit.
    tw_length : float, default=2
        Window length in seconds.
    ses_n : int, default=1
        Session number.
    plot : bool, default=True
        Whether to generate plots.
    """
    
    lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,ses_n)
    rew_ind,rew_interpolated            = ses_portpeak_DA(lick_df,photo_df,pmeta,bmeta,ses_n,tw_start,tw_length,1)
    unrew_ind,unrew_interpolated        = ses_portpeak_DA(lick_df,photo_df,pmeta,bmeta,ses_n,tw_start,tw_length,0)

    if plot:
        plot_ind_resp(rew_ind,rew_interpolated)
        plot_ind_resp(unrew_ind,unrew_interpolated)


def ses_portpeak_DA(lick_df,photo_df,pmeta,bmeta,ses_n,tw_start=0,tw_length=2,reward=1):
    """
    Extract port-specific dopamine response peaks across a session.
    
    For each port, finds the peak response to rewarded visits or the minimum
    response to unrewarded visits within a specified time window. Also returns
    event-level responses and interpolated/smoothed response traces over session
    time, sorted by port interval.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing reward status, unique-visit flags, port IDs,
        event times, and photometry indices.
    photo_df : pandas.DataFrame
        Photometry data containing a 'signal' column.
    pmeta : dict-like
        Photometry metadata containing sampling rates.
    bmeta : dict-like
        Behavioral metadata containing port intervals.
    ses_n : int
        Session number.
    tw_start : float, default=0
        Window start time relative to visit.
    tw_length : float, default=2
        Window length in seconds.
    reward : {0, 1}, default=1
        Whether to analyze rewarded or unrewarded visits.
    
    Returns
    -------
    tuple
        Event-level peak responses and interpolated/smoothed traces by port.
    """
    photo                        = photo_df['signal'].to_numpy()
    all_ind,interpolated = [],[]
    for port in range(1,7):
        photo_i    = lick_df.loc[(lick_df['rewarded']==reward) & (lick_df['unique_visit']==1) &(lick_df['port']==port)]['photo_i'].dropna()
        win_idx    = sf.extract_window(photo_i,photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)[0]
        win_signal = np.vstack([photo[win_idx[x][0]:win_idx[x][1]] for x in range(len(win_idx))])

        if reward == 1:
            peak_signal = [np.max(win_signal[i]) for i in range(len(win_signal))]

        elif reward == 0:
            peak_signal = [np.min(win_signal[i]) for i in range(len(win_signal))]
         
        # returning each individual reward response, plus event time for plotting 
        ind              = pd.DataFrame(lick_df[(lick_df['rewarded']==reward)&(lick_df['unique_visit']==1)&(lick_df['port']==port)][['event_time','photo_i']]).dropna()
        ind              = pd.DataFrame(ind['event_time'])
        ind['peak_resp'] = peak_signal
        all_ind.append(ind)
            
        # smoothing for average line 
        bin_edges  = np.arange(0,10800+1,1)
        bin_counts = np.histogram(ind['event_time'],bins=bin_edges)[0]
        # interpolate by time 
        resp_count = 0
        long_events = []
        for i,count in enumerate(bin_counts):
            if bin_counts[i]==1:
                long_events.append(ind['peak_resp'].iloc[resp_count])
                resp_count+=1
            else:
                long_events.append(np.nan)
        ind_inter = pd.DataFrame(long_events).interpolate(method='linear')
        # smooth
        ind_inter = ind_inter.rolling(window=2000,min_periods=500).mean()
        interpolated.append(ind_inter)

    all_ind         = sf.sort_data_by_interval(all_ind,bmeta['intervals'][ses_n-1])
    interpolated    = sf.sort_data_by_interval(interpolated,bmeta['intervals'][ses_n-1])

    return all_ind,interpolated 


def plot_ind_resp(all_ind,interpolated):
    """
    Plot event-level dopamine responses for each port.
    
    Creates one plot per port showing individual event responses across session
    time, with an interpolated/smoothed response trace overlaid.
    
    Parameters
    ----------
    all_ind : list of pandas.DataFrame
        Event-level responses by port. Each DataFrame should contain 'event_time'
        and 'peak_resp'.
    interpolated : list
        Smoothed response traces by port.
    """
        
    color  = plt.cm.plasma(np.linspace(0,0.9,6))
    hfont  = {'fontname':'Arial'}
    ymax   = np.ceil(np.max([all_ind[i]['peak_resp'].max() for i in range(len(all_ind))]))
    ymin   = np.floor(np.min([all_ind[i]['peak_resp'].min() for i in range(len(all_ind))]))
    for i,resps in enumerate(all_ind):
        fig,ax = plt.subplots(figsize=(2.2,1.5))
        ax.scatter(resps['event_time'],resps['peak_resp'],color=color[i],s=15,alpha=0.4,zorder=-1)
        ax.plot(range(len(interpolated[i])),interpolated[i],color='white',linewidth=4)
        ax.plot(range(len(interpolated[i])),interpolated[i],color=color[i],linewidth=2)
        ax.set_xlim([-1000,10900])
        ax.set_xticks([0,5400,10800])
        ax.set_xticklabels([0,90,180],**hfont,fontsize=18,color='black')
        ax.set_ylim([ymin-0.5,ymax])
        ax.set_yticks([])
        ax.set_yticks([0,4])
        ax.set_yticklabels([0,4],**hfont,fontsize=18,color='white')
        ax.spines[['top','right','left']].set_visible(False)
        ax.spines['bottom'].set_linewidth(1.5)
   
        ax.axhline(y=0,xmin=0.1,color='gainsboro',linestyle=(0,(1,1)),linewidth=2)
        ax.tick_params(axis='y', colors='white')

        fig.tight_layout()
