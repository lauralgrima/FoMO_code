import support_funcs as sf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pypalettes import load_cmap



def MULTIses3_DA(data_dict,tw_start=-2,tw_length=5,plot=True):
    '''
    Compute and optionally plot mean NAc dopamine responses in session 3.

    For each mouse, extracts photometry and behavioral data from session 3 and
    computes mean z-scored dopamine responses to rewarded and unrewarded visits
    for each option. Responses are aggregated across mice and can be visualized,
    with options colored by the change in expected value between sessions.

    Parameters
    ----------
    data_dict : dict
        Nested data structure containing behavioral, photometry, and metadata
        for all mice and sessions.
    tw_start : float, optional
        Start of the time window (in seconds) relative to visit/choice used for
        extracting responses. Default is -2.
    tw_length : float, optional
        Length of the time window (in seconds) for response extraction. Default is 5.
    plot : bool, optional
        If True, generates plots of mean responses across mice. Default is True.

    Returns
    -------
    None
        Computes session 3 dopamine responses and optionally generates plots.
    '''
    subset_dict = sf.subset_mice(data_dict,config=1, region='NAc')
    
    NAc_rmean_ses,NAc_umean_ses = [],[]
    
    for mouse in list(subset_dict.keys()):
        print(mouse)
        
        tot_ses   = subset_dict[mouse]['conc']['b_meta']['nsessions'] # total number of sessions for given mouse 
        ses_count = list(range(1,tot_ses+1))
        if len(ses_count)>2: # if animal did more than just the initial two sessions
            lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,3)
            rmean_ses,umean_ses = ses3_DA(lick_df,photo_df,pmeta,bmeta,3,tw_start,tw_length,first_half=False)
            if pmeta['region'] == 'NAc': # only for NAc for now 
                NAc_rmean_ses.append(rmean_ses)
                NAc_umean_ses.append(umean_ses)
                
    if plot:
        plot_ses3_DA(NAc_rmean_ses,NAc_umean_ses,tw_start,tw_length,shift=False,baseline_subtract=-1)

def ses3_DA(lick_df,photo_df,pmeta,bmeta,ses_n,tw_start=-2,tw_length=5,first_half=True):
    '''
    Compute mean dopamine responses to rewarded and unrewarded visits for each port.

    For each option, extracts photometry signal windows aligned to unique rewarded
    and unrewarded visits, then averages responses separately by reward outcome.
    Analysis can optionally be restricted to the early portion of the session.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral event data containing reward outcome, port identity, unique visit
        labels, and photometry sample indices.
    photo_df : pandas.DataFrame
        Photometry data containing the z-scored dopamine signal.
    pmeta : dict
        Photometry metadata containing sampling rate information for each session.
    bmeta : dict
        Behavioral metadata for the session.
    ses_n : int
        Session number used to index the corresponding sampling rate.
    tw_start : float, optional
        Start of the extraction window, in seconds relative to visit/choice.
        Default is -2.
    tw_length : float, optional
        Length of the extraction window, in seconds. Default is 5.
    first_half : bool, optional
        If True, limits averaging to the early portion of the session. Default is True.

    Returns
    -------
    rmean_ses : pandas.DataFrame
        Mean photometry responses to rewarded visits for each port.
        Rows correspond to ports and columns correspond to time samples.
    umean_ses : pandas.DataFrame
        Mean photometry responses to unrewarded visits for each port.
        Rows correspond to ports and columns correspond to time samples.
    '''
    
    photo = photo_df['signal'].to_numpy()
    # session 3 - rewarded trials 
    rmean_ses,umean_ses = [],[]
    for port in range(1,7):
        rphoto_i    = lick_df.loc[(lick_df['rewarded']==1) & (lick_df['unique_visit']==1)&(lick_df['port']==port)]['photo_i'].dropna()
        rwin_idx    = sf.extract_window(rphoto_i,photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)[0]
        rwin_signal = np.vstack([photo[rwin_idx[x][0]:rwin_idx[x][1]] for x in range(len(rwin_idx))])
        if first_half:
            rwin_signal = rwin_signal[:np.round(len(rwin_signal)/4).astype(int)]
        rmean_ses.append(np.mean(rwin_signal,axis=0))
        
        uphoto_i    = lick_df.loc[(lick_df['rewarded']==0) & (lick_df['unique_visit']==1)&(lick_df['port']==port)]['photo_i'].dropna()
        uwin_idx    = sf.extract_window(uphoto_i,photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)[0]
        uwin_signal = np.vstack([photo[uwin_idx[x][0]:uwin_idx[x][1]] for x in range(len(uwin_idx))])
        if first_half:
            uwin_signal = uwin_signal[:np.round(len(uwin_signal)/4).astype(int)]
        umean_ses.append(np.mean(uwin_signal,axis=0))
        
    # if a certain trial type doesn't occur for a port 
    for i,sig in enumerate(rmean_ses):
        if rmean_ses[i].size == 0:
            rmean_ses[i] = np.array([np.nan]*len(rmean_ses[i]))
        if umean_ses[i].size == 0:
            umean_ses[i] = np.array([np.nan]*len(umean_ses[i]))
            
    return pd.DataFrame(rmean_ses),pd.DataFrame(umean_ses)
                 
   
def plot_ses3_DA(rew_mean,unrew_mean,tw_start,tw_length,shift=0.3,baseline_subtract=-1):
    '''
    Plot mean session 3 dopamine responses by reward outcome and option.

    Generates separate plots for rewarded and unrewarded visits. For each outcome,
    responses are averaged across mice and plotted for each option. Traces are
    colored by the change in expected value from session 2 to session 3, with
    greener colors indicating larger positive changes in expectation.

    Parameters
    ----------
    rew_mean : list of pandas.DataFrame
        Mean rewarded-visit responses for each mouse. Rows correspond to ports and
        columns correspond to time samples.
    unrew_mean : list of pandas.DataFrame
        Mean unrewarded-visit responses for each mouse. Same structure as `rew_mean`.
    tw_start : float
        Start of the plotted time window, in seconds relative to visit/choice.
    tw_length : float
        Length of the plotted time window, in seconds.
    shift : float or bool, optional
        Vertical offset applied between port traces for visualization. If False,
        traces are plotted without vertical offset. Default is 0.3.
    baseline_subtract : float or bool, optional
        Time point, in seconds relative to visit/choice, used for baseline subtraction.
        If False, no baseline subtraction is applied. Default is -1.

    Returns
    -------
    None
        Generates matplotlib figures.
    '''
    # averaging
    Rses3_mean = pd.concat(rew_mean,axis=0).groupby(pd.concat(rew_mean,axis=0).index).mean()
    Uses3_mean = pd.concat(unrew_mean,axis=0).groupby(pd.concat(unrew_mean,axis=0).index).mean()
    means      = [Rses3_mean,Uses3_mean]
    
    # plotting
    hfont = {'fontname':'Arial'} 
    
    x       = list(range(0,len(Rses3_mean.iloc[0])))
    x_sec   = [sample/130 for sample in x] # turn samples into seconds 
    x_old   = [0,0-tw_start,tw_length]
    x_ticks = [tw_start,0,tw_length--tw_start]
    
    # colors - based on change in expectation (positive or negative)
    cmap      = load_cmap("ArmyRose")
    diff_exp  = np.array([-5,-1.5,-1.5,2.5,1.5,4])
    scale_exp = (diff_exp-diff_exp.min())/(diff_exp.max()-diff_exp.min())
    colors    = [cmap(i) for i in scale_exp][::-1]

    for i in [0,1]: # for rewarded and unrewarded trials 
        fig,ax = plt.subplots(figsize=(3,3))
        
        for port in range(1,7):
            port_sig = means[i].iloc[port-1]
            if baseline_subtract:
                baseline_samp = abs(tw_start)*130-abs(baseline_subtract)*130
                port_sig      = port_sig - port_sig[baseline_samp]
            if shift:
                port_sig = port_sig-(port*shift)
            ax.plot(x_sec,port_sig,color=colors[port-1],label=port)
            ax.set_xlabel('time from visit (s)',**hfont,fontsize=18)
            if shift: 
                ax.axhline(y=0-(port*shift),color='dimgrey',linewidth=1.5,alpha=0.3,xmin=0.05,xmax=0.95)
                
        ax.axvline(x=abs(tw_start),color='black',linestyle='--',linewidth=1.5)
        ax.set_xticks(x_old)
        ax.set_xticklabels(x_ticks)
        if i == 0:
            ax.set_yticks([0,1])
        else:
            ax.set_yticks([0,0.3])
        ax.spines[['right','top','left']].set_visible(False)
        ax.spines[['bottom']].set_position(('outward', 10))
        ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
        ax.yaxis.set_tick_params(width=1.5)
        ax.xaxis.set_tick_params(width=1.5)
        ax.spines[['left','bottom']].set_linewidth(1.5)
        
        fig.tight_layout()
