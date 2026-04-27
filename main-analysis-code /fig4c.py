import support_funcs as sf
import numpy as np
import matplotlib.pyplot as plt


def MULTIvisit_DA_by_port(data_dict,ses_n=1,tw_start=-2,tw_length=5,plot=True):
    """
    Plot visit-aligned dopamine responses by port.
    
    Plots one average response trace per port, with colors indicating port identity.
    
    Parameters
    ----------
    tw_start : float
        Window start time relative to visit.
    tw_length : float
        Window length in seconds.
    signal : list
        List of response traces, one per port.
    """
    for region in ['NAc','DMS']: # plotting each region separately 
        subset_dict = sf.subset_mice(data_dict,config=1,region=region)
        
        all_rew,all_unrew = [],[]
        for mouse in subset_dict:
            print(mouse)

            lick_df,ses_photo_df,vid_df,bmeta,pmeta   = sf.extract_data(data_dict,mouse,ses_n)
            win_signal_rew_port,win_signal_unrew_port = port_DA(lick_df,ses_photo_df,pmeta,bmeta,ses_n,tw_start=tw_start,tw_length=tw_length)

            all_rew.append(win_signal_rew_port)
            all_unrew.append(win_signal_unrew_port)
            
        # average for each port 
        av_rew_port   = np.mean(np.array(all_rew),axis=0)
        av_unrew_port = np.mean(np.array(all_unrew),axis=0)

        if plot:
            plot_rew_unrew_port(tw_start,tw_length,av_rew_port)
            plot_rew_unrew_port(tw_start,tw_length,av_unrew_port)
            
            
            

def port_DA(lick_df,photo_df,pmeta,bmeta,ses_n,tw_start=-2,tw_length=5):
    """
    Extract average dopamine responses by port.
    
    For each port, extracts visit-aligned photometry windows separately for
    rewarded and unrewarded unique visits, averages responses within port, and
    sorts ports by interval.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing reward status, unique-visit flags, port IDs,
        and photometry indices.
    photo_df : pandas.DataFrame
        Photometry data containing a 'signal' column.
    pmeta : dict-like
        Photometry metadata containing sampling rates.
    bmeta : dict-like
        Behavioral metadata containing port intervals.
    ses_n : int
        Session number.
    tw_start : float, default=-2
        Window start time relative to visit.
    tw_length : float, default=5
        Window length in seconds.
    
    Returns
    -------
    tuple
        Rewarded and unrewarded mean response traces by port, sorted by interval.
    """
    photo = photo_df['signal'].to_numpy()
    
    rew_by_port,unrew_by_port = [],[]
    for port in range(1,7):
        photo_i_rew        = lick_df.loc[(lick_df['rewarded']==1) & (lick_df['unique_visit']==1)][['photo_i','port']].dropna()
        photo_i_unrew      = lick_df.loc[(lick_df['rewarded']==0) & (lick_df['unique_visit']==1)][['photo_i','port']].dropna()
        photo_i_rew_port   = photo_i_rew[photo_i_rew['port']==port]['photo_i']
        photo_i_unrew_port = photo_i_unrew[photo_i_unrew['port']==port]['photo_i']
            
        if len(photo_i_rew_port)>0:
                
            win_idx_rew_port      = sf.extract_window(photo_i_rew_port,photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)[0]
            win_idx_unrew_port    = sf.extract_window(photo_i_unrew_port,photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)[0]
            
            win_signal_rew_port   = np.vstack([photo[win_idx_rew_port[x][0]:win_idx_rew_port[x][1]] for x in range(len(win_idx_rew_port))])
            win_signal_unrew_port = np.vstack([photo[win_idx_unrew_port[x][0]:win_idx_unrew_port[x][1]] for x in range(len(win_idx_unrew_port))])
            
            av_win_sig_rew_port   = np.mean(win_signal_rew_port,axis=0)
            av_win_sig_unrew_port = np.mean(win_signal_unrew_port,axis=0)
                
            rew_by_port.append(av_win_sig_rew_port)
            unrew_by_port.append(av_win_sig_unrew_port)
    
    # sort by interval
    rew_by_port   = sf.sort_data_by_interval(rew_by_port,np.array(bmeta['intervals'][ses_n-1]))
    unrew_by_port = sf.sort_data_by_interval(unrew_by_port,np.array(bmeta['intervals'][ses_n-1]))

    return rew_by_port,unrew_by_port


def plot_rew_unrew_port(tw_start,tw_length,signal):    
    """
    Plot visit-aligned dopamine responses for each port.
    
    Generates a line plot with one trace per port (1–6), showing the average
    photometry signal aligned to visit time.
    
    Parameters
    ----------
    tw_start : float
        Window start time relative to visit.
    tw_length : float
        Window length in seconds.
    signal : list
        List of length 6, where each element is a 1D array representing the
        average response trace for a port.
    """
    hfont = {'fontname':'Arial'} 
    fig,ax = plt.subplots(figsize=(3.5,3.5))

    colors = plt.cm.plasma(np.linspace(0,0.9,6))

    for port in range(1,7):

        port_sig = signal[port-1]
        
        # averaging 
        x       = list(range(0,len(port_sig)))
        x_sec   = [sample/130 for sample in x] # turn samples into seconds 
        x_old   = [0,0-tw_start,tw_length]
        x_ticks = [tw_start,0,tw_length--tw_start]
    
        ax.plot(x_sec,port_sig,color=colors[port-1],label=port)
        ax.set_xlabel('time from visit (s)',**hfont,fontsize=18)
        ax.axhline(y=0,color='gainsboro',linewidth=2,alpha=0.3,xmin=0.05,xmax=0.95,linestyle=(0,(1,1)))
      
    ax.axvline(x=abs(tw_start),color='gainsboro',linewidth=2, zorder=-1)
    ax.set_xticks(x_old)
    ax.set_xticklabels(x_ticks)
    ax.set_ylim([-1,2])
    ax.set_yticks([0,0.5])
    ax.spines[['right','top','left']].set_visible(False)
    ax.spines[['bottom']].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.yaxis.set_tick_params(width=1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.spines[['left','bottom']].set_linewidth(1.5)

    fig.tight_layout()
