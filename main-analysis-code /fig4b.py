import support_funcs as sf
import numpy as np
import matplotlib.pyplot as plt
import stats

def MULTIvisit_DA(data_dict,ses_n=[1,2],tw_start=-2,tw_length=5,plot=True):
    """
    Plot average dopamine responses to rewarded and unrewarded visits across mice.
    
    For each region, extracts visit-aligned photometry responses across selected
    sessions, averages within mouse, and optionally plots rewarded and unrewarded
    response traces.
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral and photometry data.
    ses_n : list, default=[1, 2]
        Sessions to include.
    tw_start : float, default=-2
        Window start time relative to visit.
    tw_length : float, default=5
        Window length in seconds.
    plot : bool, default=True
        Whether to plot average traces for each region.
    """

    for region in ['NAc','DMS']: # plotting each region separately 
        subset_dict = sf.subset_mice(data_dict,config=1,region=region)
        
        all_rew,all_unrew = [],[]
        for mouse in subset_dict:
            print(mouse)
            
            tot_ses   = data_dict[mouse]['conc']['b_meta']['nsessions'] # total number of sessions for given mouse 
            ses_count = list(range(1,tot_ses+1))
            sessions  = list(set(ses_count).intersection(ses_n))
            if sessions: 
                arew,aunrew = [],[]
                for ses in sessions:
                    lick_df,ses_photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,ses)
                    win_signal_rew,win_signal_unrew         = visit_DA(lick_df,ses_photo_df,pmeta,ses,tw_start=tw_start,tw_length=tw_length)
                    arew.append(win_signal_rew), aunrew.append(win_signal_unrew)
                
                mrew   = np.mean(arew,axis=0)
                munrew = np.mean(aunrew,axis=0)
                
            all_rew.append(mrew)
            all_unrew.append(munrew)
            
        if plot:
            plot_rew_unrew_resp(tw_start,tw_length,[all_rew,all_unrew])
            
            
def visit_DA(lick_df,photo_df,pmeta,ses_n,tw_start=-2,tw_length=5):
    """
    Extract average dopamine signal around rewarded and unrewarded visits.
    
    For unique visits, extracts photometry signal windows around visit-aligned
    photometry indices and returns the mean traces separately for rewarded and
    unrewarded visits.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing reward status, unique-visit flags, port IDs,
        and photometry indices.
    photo_df : pandas.DataFrame
        Photometry data containing a 'signal' column.
    pmeta : dict-like
        Photometry metadata containing sampling rates.
    ses_n : int
        Session number.
    tw_start : float, default=-2
        Window start time relative to event.
    tw_length : float, default=5
        Window length.
    
    Returns
    -------
    tuple of numpy.ndarray
        Mean rewarded-visit trace and mean unrewarded-visit trace.
    """
    
    photo = photo_df['signal'].to_numpy()
    
    photo_i_rew   = lick_df.loc[(lick_df['rewarded']==1) & (lick_df['unique_visit']==1)][['photo_i','port']].dropna()
    photo_i_unrew = lick_df.loc[(lick_df['rewarded']==0) & (lick_df['unique_visit']==1)][['photo_i','port']].dropna()

    win_idx_rew      = sf.extract_window(photo_i_rew['photo_i'],photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)[0]
    win_idx_unrew    = sf.extract_window(photo_i_unrew['photo_i'],photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)[0]

    win_signal_rew      = np.vstack([photo[win_idx_rew[x][0]:win_idx_rew[x][1]] for x in range(len(win_idx_rew))])
    win_signal_unrew    = np.vstack([photo[win_idx_unrew[x][0]:win_idx_unrew[x][1]] for x in range(len(win_idx_unrew))])

    return np.mean(win_signal_rew,axis=0),np.mean(win_signal_unrew,axis=0)


            
### PLOTTING                    
            
def plot_rew_unrew_resp(tw_start,tw_length,all_signal):
    """
    Plot average photometry responses for rewarded and unrewarded visits.
    
    Plots mean traces with SEM shading for two visit conditions, aligned to visit
    time.
    
    Parameters
    ----------
    tw_start : float
        Window start time relative to visit.
    tw_length : float
        Window length in seconds.
    all_signal : list
        List containing rewarded and unrewarded signal arrays across sessions/mice.
    """
    hfont = {'fontname':'Helvetica'} 
    fig,ax = plt.subplots(1,1,figsize=(3.5,3.5))

    x       = list(range(0,len(all_signal[0][0])))
    x_sec   = [sample/130 for sample in x] # turn samples into seconds 
    x_old   = [0,0-tw_start,tw_length]
    x_ticks = [tw_start,0,tw_length--tw_start]

    # averaging
    maxes    = []
    pl_lines = []
    color    = ['navy','firebrick']
    for i,signal in enumerate(all_signal):
        sig = [sig for sig in all_signal[i] if not np.isnan(sig).any()]
        ses_signal_mean = np.nanmean(sig,axis=0)
        ses_signal_sem  = stats.sem(sig,axis=0,nan_policy='omit')
        ses_pos_sem     = ses_signal_mean+ses_signal_sem
        ses_neg_sem     = ses_signal_mean-ses_signal_sem
        
        maxes.append(np.max(ses_signal_mean))

        pl_line = ax.plot(x_sec,ses_signal_mean,color=color[i])
        ax.fill_between(x_sec,ses_pos_sem,ses_neg_sem,alpha=0.3,color=color[i])
        pl_lines.append(pl_line)

    ax.axvline(x=abs(tw_start),color='gainsboro',linewidth=2,zorder=-1)
    ax.axhline(y=0,color='gainsboro',linestyle=(0,(1,1)),linewidth=2,zorder=-1)
    ax.set_xticks(x_old)
    ax.set_xticklabels(x_ticks)
    ax.set_ylim([-1,2])
    ax.set_yticks([0,1])
    ax.set_xlim([0,5])
    ax.spines[['right','top','left']].set_visible(False)
    ax.spines[['bottom']].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.yaxis.set_tick_params(width=1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    
    ax.set_xlabel('time from visit (s)',**hfont,fontsize=18)

    fig.tight_layout()
    plt.show()


