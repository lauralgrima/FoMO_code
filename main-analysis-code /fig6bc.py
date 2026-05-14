import os 
import numpy as np
import support_funcs as sf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import import_mat as im
from pypalettes import load_cmap

def MULTIses_pchoice(data_dict,mat_GLM_path,win_len=20,plot=True,plot_example=True):
    '''
    Compute choice probabilities across sessions for behavior and AQUA model data.

    For each mouse in the NAc and DMS config 1 cohorts, extracts behavioral visits
    and corresponding AQUA model visits across sessions, computes the fraction of
    choices to each port in non-overlapping time windows, and concatenates results
    across sessions. NAc mice are additionally stored separately for example plotting.

    Parameters
    ----------
    data_dict : dict
        Nested behavioral/photometry data dictionary containing mouse and session data.
    mat_GLM_path : str
        Path to MATLAB GLM/AQUA model output files.
    win_len : int or float, optional
        Window length in minutes used to compute choice probabilities. Default is 20.
    plot : bool, optional
        If True, plots across-session choice probability trajectories. Default is True.

    Returns
    -------
    None
        Computes choice probability data and optionally generates plots.
    '''

    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    subset_dict     = {**subset_dict_NAc,**subset_dict_DMS}
        
    _,AQUA_beh_NAc  = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS  = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_dicts      = {**AQUA_beh_DMS,**AQUA_beh_NAc}
    
    all_props,all_AQUA,NAc_props,NAc_AQUA_props,NAc_mice = [],[],[],[],[]
    for mouse in list(subset_dict.keys()):
        anm_props,AQUA_props = [],[]
        print(mouse)
        
        tot_ses   = subset_dict[mouse]['conc']['b_meta']['nsessions'] # total number of sessions for given mouse 
        ses_count = list(range(1,tot_ses+1))
        
        # extract relevant model data 
        mAQUA = AQUA_dicts[mouse]
    
        if len(ses_count)>2: # if animal did more than just the initial two sessions
            for ses in ses_count:
                print(ses)
                lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,ses)
                # get visits
                beh_visits                          = lick_df[lick_df['unique_visit']==1][['port','event_time']]
                AQUA_visits                         = im.get_visits_AQUA(mAQUA['AQUA_vis_invis'],lick_df,ses,bmeta)   
                # call prop function
                anm_props.append(pchoice(beh_visits,bmeta,win_len=win_len))
                AQUA_props.append(pchoice(AQUA_visits,bmeta,win_len=win_len))
            all_anm_props  = pd.concat(anm_props).reset_index(drop=True)
            all_AQUA_props = pd.concat(AQUA_props).reset_index(drop=True)
            all_props.append(all_anm_props)
            all_AQUA.append(all_AQUA_props)
            
            if pmeta['region'] == 'NAc':
                NAc_mice.append(mouse)
                NAc_props.append(all_anm_props)
                NAc_AQUA_props.append(all_AQUA_props)

        elif len(ses_count)<2:
            continue 
        
    if plot: 
        # fig 6b
        plot_props(all_props,win_len,smooth=False)
        plot_props(all_AQUA,win_len,smooth=False)
        
    if plot_example:
        # fig 6c
        plot_example_prop(NAc_props,NAc_AQUA_props,NAc_mice)
        
        

def pchoice(visits,bmeta,win_len):
    '''
    Calculate the fraction of choices made to each port over time.

    Choice fractions are computed in non-overlapping time windows across a
    3-hour session. Values are returned in port order, with port columns
    reordered for the 1200 s interval configuration to match the standard
    plotting/order convention.

    Parameters
    ----------
    visits : pandas.DataFrame
        Visit-level behavioral data. Must contain 'event_time' and 'port' columns.
    bmeta : dict
        Behavioral metadata containing interval configuration information.
    win_len : int or float
        Window length in minutes.

    Returns
    -------
    props : pandas.DataFrame
        Fraction of choices to each port within each time window.
        Rows correspond to time windows and columns correspond to ports.
    '''
    win_lenm = win_len*60
    windows  = np.arange(0,10800+win_lenm,win_lenm)
    win_vis  = [visits[(visits['event_time']>windows[i])&(visits['event_time']<windows[i+1])] for i in range(len(windows)-1)]
    
    props    = [np.array(win_data['port'].value_counts().sort_index()/len(win_data)) for win_data in win_vis] # in port order (not port rank)
    props    = pd.DataFrame(props).fillna(0)
    
    if bmeta['intervals'][0][0] == 1200:
        props = props[[4,5,3,2,0,1]]
    elif bmeta['intervals'][0][0] == 30:
        props = props
    else:
        print('other config')

    return props


### PLOTTING

def plot_example_prop(NAc_props,NAc_AQUA_props,NAc_mice):
    """
    Plot example choice probability trajectories for individual mice.
    
    For each mouse, plots the probability of sampling each port (0–5) across time,
    with ports overlaid in different colors. Behavior and corresponding model data
    are plotted separately. Vertical dashed lines indicate session boundaries.
    
    Parameters
    ----------
    NAc_props : list of pandas.DataFrame
        Behavioral choice probability data for each mouse (columns = ports).
    NAc_AQUA_props : list of pandas.DataFrame
        Model-derived choice probabilities for each mouse (same structure as NAc_props).
    NAc_mice : list of str
        Mouse identifiers used for plot titles and file naming.
    savefig : str or bool
        If a string, directory path to save figures. If False, figures are not saved.
    
    Returns
    -------
    None
        Generates and optionally saves matplotlib figures.
    """
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    port_cols = [load_cmap("PassionFruit")(i) for i in range(0,6)]
    ses_cols  = [load_cmap("Starfish")(i) for i in [0.25,0.5,0.75,1.0]]

    for j,datatype in enumerate([NAc_props,NAc_AQUA_props]):
        for i,pmouse in enumerate(datatype):
            # plotting mouse data
            fig,ax = plt.subplots(figsize=(3,3.5))
            for port in range(0,6):
                option = pmouse.columns.get_loc(port)
                ax.plot(pmouse[port].rolling(window=5).mean().dropna(),color=port_cols[option])
                
                ax.axvline(x=9, ymin=0,ymax=1,linestyle='dotted',color=ses_cols[0])
                ax.axvline(x=18,ymin=0,ymax=1,linestyle='dotted',color=ses_cols[1])
                ax.axvline(x=27,ymin=0,ymax=1,linestyle='dotted',color=ses_cols[2])
                ax.axvline(x=36,ymin=0,ymax=1,linestyle='dotted',color=ses_cols[3])
                
                ax.spines[['left','top','right','bottom']].set_visible(False)
                ax.xaxis.set_ticks([])
                ax.xaxis.set_tick_params(length=0)
                if j == 0:
                    ax.set_title(NAc_mice[i]+'beh')
                elif j == 1:
                    ax.set_title(NAc_mice[i]+'model')
                ax.yaxis.set_ticks([])
                ax.yaxis.set_tick_params(length=0)
                
                ax.set_ylim([0,0.4])
                
                fig.tight_layout() 


def plot_props(all_props,win_len,smooth):
    """
    Plot port choice probabilities across time alongside changing port reward qualities.
    
    For each of six ports, plots individual trajectories (gray) and the across-subject median (black).
    Port reward quality is overlaid on a secondary axis and segmented by session.
    
    Parameters
    ----------
    all_props : list of pandas.DataFrame
        Per-subject dataframes containing choice probabilities (columns = ports).
    win_len : int or float
        Window length in minutes (used for defining session structure).
    smooth : bool
        If True, applies a rolling mean (window=2) to smooth trajectories.
    """
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    
    if smooth:
        aprops = [all_props[i].rolling(window=2).mean().dropna().reset_index(drop=True) for i in range(len(all_props))]
    else:
        aprops = all_props

    port_means = [np.array(pd.concat([aprops[i].iloc[:,j] for i in range(len(aprops))],axis=1).median(axis=1)) for j in range(0,6)]
    
    n_x = len(port_means[0])
    x_vals = np.arange(n_x)

    pq = [[1,1,6,3.5,2],
          [2,2,3.5,5,5],
          [3.5,3.5,5,2,1],
          [3.5,3.5,1,1,6],
          [5,5,3.5,3.5,3.5],
          [6,6,2,6,3.5]]

    pq_line_col = plt.cm.plasma(np.linspace(0,0.9,6))
    
    for port in range(0,6):
        fig,ax = plt.subplots(figsize=(2,3))
        ax2    = ax.twinx()
        ax2.invert_yaxis()

        # plot individual traces
        for props in aprops:
            ax.plot(np.arange(len(props)),props.iloc[:,port].values,color='gainsboro',linewidth=1.5)
            
        # plot median
        ax.plot(x_vals,port_means[port],color='black',linewidth=1.5)

        # plot port quality across same x axis
        session_edges = np.linspace(0,n_x-1,len(pq[port])+1)

        for i,qual in enumerate(pq[port]):
            icol = qual-1
            if icol == 2.5:
                icol = 3
            icol = int(icol)

            ax2.plot([session_edges[i],session_edges[i+1]],[qual,qual],color=pq_line_col[icol])

        ax.set_xlim([0,n_x-1])
        ax.set_ylim([0,0.5])
        ax2.set_ylim([6.5,0.5])
        
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=True)
        ax.set_xticks([])

        ax.spines[['left','top','right','bottom']].set_visible(False)
        ax2.spines[['left','top','right','bottom']].set_visible(False)
        
        cmap    = load_cmap("Starfish")
        sescolors = [cmap(i) for i in [0.25,0.5,0.75,1.0]]

        for x,c in zip(session_edges[1:-1],sescolors):
            ax.axvline(x=x,ymin=0,ymax=1,linestyle='dotted',color=c)
        
        fig.tight_layout()

