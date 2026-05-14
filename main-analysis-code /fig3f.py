import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import support_funcs as sf
import import_mat as im
from scipy.stats import entropy

def MULTI_kl_withinses(data_dict,mat_GLM_path,ses_n=1,win_len=50,max_vis=300,ses_fraction=4,portion='end',plot=True):
    """
    Compute within-session KL divergence for mouse behavior and AQUA model data.
    
    For each mouse, estimates a baseline port distribution from a specified session
    portion, then computes sliding-window KL divergence across the session. Results
    are optionally plotted for both experimental and AQUA-derived behavior.
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral data.
    mat_GLM_path : str
        Path to GLM/AQUA model data.
    ses_n : int, default=1
        Session number to analyze.
    win_len : int, default=50
        Sliding window size for KL divergence.
    max_vis : int, default=300
        Maximum number of visits to include.
    ses_fraction : int, default=4
        Fraction denominator for baseline period.
    portion : {'end', 'beginning'}, default='end'
        Session portion used to estimate baseline probabilities.
    plot : bool, default=True
        Whether to plot KL divergence traces.
    """

    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    subset_dict     = {**subset_dict_NAc,**subset_dict_DMS}
    
    _,AQUA_beh_NAc  = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS  = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_beh        = {**AQUA_beh_DMS,**AQUA_beh_NAc}
    
    mouse_kl,AQUA_kl = [],[]
    
    for mouse in subset_dict:
        # mouse data
        print(mouse)
        lick_df,_,_,bmeta,_ = sf.extract_data(data_dict,mouse,ses_n)  
        port_probs          = baseline_prob_dist(lick_df,ses_fraction,data_type='real',portion=portion) # port probs at baseline
        kl_div              = kl_divergence_within(lick_df,port_probs,win_len=win_len,data_type='real')
        mouse_kl.append(pd.DataFrame(kl_div).rolling(window=win_len,min_periods=2).mean())
        
        # model data
        AQUA_visits           = im.get_visits_AQUA(AQUA_beh[mouse]['AQUA_vis_invis'],lick_df,ses_n,bmeta)['port']
        AQUA_visits           = AQUA_visits.map(sf.port_2_rank(bmeta,ses_n))
        AQUA_port_probs       = baseline_prob_dist(AQUA_visits,ses_fraction,data_type='AQUA',portion=portion) # port probs at baseline
        kl_div                = kl_divergence_within(AQUA_visits,AQUA_port_probs,win_len=win_len,data_type='AQUA')
        AQUA_kl.append(pd.DataFrame(kl_div).rolling(window=win_len,min_periods=2).mean())
        
    exp_kl      = pd.concat(mouse_kl,axis=1).iloc[1:max_vis]
    exp_AQUA_kl = pd.concat(AQUA_kl,axis=1).iloc[1:max_vis]

    if plot: # plot one session at a time
        plot_within_kl(exp_kl,max_vis)
        plot_within_kl(exp_AQUA_kl,max_vis)
            

def kl_divergence_within(lick_df,port_probs,win_len=50,data_type='real'):
    """
    Compute the Kullback–Leibler (KL) divergence over a moving window within a session.
    
    This function compares the empirical distribution of port visits within sliding
    windows to a baseline distribution (`port_probs`). The KL divergence is computed
    for each window, producing a time series of divergence values across the session.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame or pandas.Series
        Input data containing port visit information. If `data_type='real'`, this
        must be a DataFrame with columns:
            - 'unique_visit': indicator (1 for first lick in a visit)
            - 'port': port identity
        If `data_type!='real'`, `lick_df` is treated directly as a sequence of ports.
    
    port_probs : array-like
        Baseline probability distribution over ports. Must be aligned with the
        output of `conv_to_prop`.
    
    win_len : int, optional (default=50)
        Size of the sliding window used to compute local port distributions.
    
    data_type : str, optional (default='real')
        Determines how input data is interpreted:
            - 'real': filter `lick_df` to include only unique visits
            - otherwise: use `lick_df` directly as a sequence
    
    Returns
    -------
    numpy.ndarray
        Array of KL divergence values, one per valid window, representing how much
        each local window distribution deviates from the baseline distribution.
    
    Notes
    -----
    - Uses a sliding window approach to compute local port probabilities via
      `conv_to_prop`.
    - KL divergence is computed using `scipy.stats.entropy`.
    - Windows that cannot form a valid probability distribution (e.g., insufficient
      data) are dropped.
    """
    if data_type == 'real':
        visitsp = lick_df[lick_df['unique_visit']==1]['port']
    else:
        visitsp = lick_df 

    port_props_mov = [conv_to_prop(visitsp.iloc[i:win_len+i]) for i in range(len(visitsp))]
    mov_props      = pd.concat(port_props_mov,axis=1).transpose().dropna()
    kl_from_ses    = np.array([entropy(port_probs,mov_props.iloc[i,:]) for i in range(0,len(mov_props))])

    return kl_from_ses


def baseline_prob_dist(lick_df,ses_fraction,data_type='real',portion='end'):
    """
    Estimate a baseline port probability distribution from a session subset.
    
    Selects a fraction of visits from either the beginning or end of the session
    and computes the normalized frequency of port choices.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame or sequence
        If `data_type='real'`, expects columns 'unique_visit' and 'port' and uses
        only unique visits. Otherwise treated as a sequence of ports.
    ses_fraction : int
        Defines the fraction of the session to use (e.g., 4 = last/first quarter).
    data_type : str, default='real'
        Controls how `lick_df` is interpreted.
    portion : {'end', 'beginning'}, default='end'
        Whether to use the last or first fraction of the session.
    
    Returns
    -------
    pandas.Series
        Probability distribution over ports. Missing ports are assigned a small
        nonzero probability.
    """
    if data_type == 'real':
        visits = lick_df[lick_df['unique_visit']==1]['port']
    else:
        visits = lick_df 
        
    if portion == 'end':
        vis_fract   = visits.iloc[len(visits)-np.round(len(visits)/ses_fraction).astype(int):]
    elif portion == 'beginning':
        vis_fract   = visits.iloc[:len(visits)-(len(visits)-np.round(len(visits)/ses_fraction)).astype(int)]
    port_probs  = vis_fract.value_counts().sort_index()/np.sum(vis_fract.value_counts().sort_index())
    if len(port_probs)<6: # option wasn't sampled in that timeframe 
        missing_opt = list(set(np.array([1,2,3,4,5,6]))-set(np.array(port_probs.index)))
        port_probs  = pd.concat([port_probs,pd.DataFrame([0.0001],index=missing_opt)]).sort_index()
        port_probs  = pd.Series(port_probs.iloc[:,0])
    
    return port_probs 



def conv_to_prop(sub_visits): 
    """
    Convert a sequence of port visits into a probability distribution.
    
    Computes normalized frequencies of visits to each port, ensuring that all
    possible ports (1–6) are represented. Ports not visited in the input are
    assigned zero probability.
    
    Parameters
    ----------
    sub_visits : array-like
        Sequence of port visits.
    
    Returns
    -------
    pandas.Series
        Probability distribution over ports (indexed 1–6).
    """
    sub_visits  = pd.Series(sub_visits)
    missing_opt = np.array(list(set(np.array([1,2,3,4,5,6]))-set(sub_visits.unique())))
    vis_props   = sub_visits.value_counts().sort_index()/np.sum(sub_visits.value_counts())
    if len(missing_opt)>0:
        for miss in missing_opt:
            vis_props = pd.concat([vis_props,pd.DataFrame([0],index=[miss])]).sort_index()
    return vis_props            
            

### PLOTTING

def plot_within_kl(exp_kl,max_vis):
    """
    Plot within-session KL divergence across mice.
    
    Each individual trace represents one mouse, and the bold black trace shows
    the across-mouse average.
    
    Parameters
    ----------
    exp_kl : pandas.DataFrame
        KL divergence values, with visits along rows and mice along columns.
    max_vis : int
        Maximum visit number to display on the x-axis.
    """
    mean_kl = pd.DataFrame(exp_kl).mean(axis=1)

    hfont   = {'fontname':'Arial'} 
    fig, ax = plt.subplots(figsize=(3.2,2.5))

    for i in range(len(exp_kl.iloc[0])):
        ax.plot(exp_kl.iloc[:,i],color='gainsboro',linewidth=1.5)
    ax.plot(mean_kl,color='black',linewidth=1.5)

    ax.set_ylim([0,0.6])
    ax.set_yticks([0,0.3,0.6])
    ax.set_yticklabels([0,0.3,0.6],**hfont,fontsize=18)

    ax.set_xlim([0,max_vis])
    ax.set_xticks([0,max_vis])
    ax.set_xticklabels([0,max_vis],**hfont,fontsize=18)
    ax.set_xlabel('visit number', **hfont,fontsize=18)
    ax.set_ylabel('KL divergence (bits)', **hfont,fontsize=18)
    ax.spines[['right','top']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')

    plt.show()
    fig.tight_layout()





