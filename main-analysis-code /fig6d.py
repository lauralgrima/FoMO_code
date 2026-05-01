import pandas as pd
import numpy as np
import support_funcs as sf 
import matplotlib.pyplot as plt
import import_mat as im
from scipy.stats import entropy

def MULTImultises_kl(data_dict,mat_GLM_path,lim2=200,lim3=300,win_len=50,plot=True):
    '''
    Calculate KL divergence around the session 2 to session 3 transition.

    For each mouse, compares port choice distributions around the session switch
    to a baseline distribution estimated from the final quarter of session 2.
    KL divergence is calculated for the last `lim2` choices of session 2 and the
    first `lim3` choices of session 3 using sliding windows of length `win_len`.
    Behavioral data and corresponding AQUA model data are analyzed separately.

    Parameters
    ----------
    data_dict : dict
        Nested behavioral/photometry data dictionary containing mouse and session data.
    mat_GLM_path : str
        Path to MATLAB GLM/AQUA model output files.
    lim2 : int, optional
        Number of choices from the end of session 2 to include. Default is 200.
    lim3 : int, optional
        Number of choices from the start of session 3 to include. Default is 300.
    win_len : int, optional
        Sliding window length, in choices, used to estimate choice distributions
        for KL divergence. Default is 50.
    plot : bool, optional
        If True, plots KL divergence traces for behavioral and model data.

    Returns
    -------
    None
        Computes KL divergence traces and optionally generates plots.
    '''
    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    subset_dict     = {**subset_dict_NAc,**subset_dict_DMS}
        
    _,AQUA_beh_NAc  = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS  = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_dicts      = {**AQUA_beh_DMS,**AQUA_beh_NAc}
    
    across_ses_kl,across_ses_kl_model = [],[]
    for mouse in list(subset_dict.keys()):
        print(mouse)
        tot_ses   = data_dict[mouse]['conc']['b_meta']['nsessions'] # total number of sessions for given mouse 
        ses_count = list(range(1,tot_ses+1))
        mAQUA     = AQUA_dicts[mouse]
        
        if len(ses_count)>2:
            # session 2
            lick_df2,_,_,bmeta,_ = sf.extract_data(subset_dict,mouse,2) 
            AQUA_visits2         = im.get_visits_AQUA(mAQUA['AQUA_vis_invis'],lick_df2,2,bmeta)
            # session 3
            lick_df3,_,_,bmeta,_ = sf.extract_data(subset_dict,mouse,3) 
            AQUA_visits3         = im.get_visits_AQUA(mAQUA['AQUA_vis_invis'],lick_df3,3,bmeta)
            
            ses2_pps             = baseline_prob_dist(lick_df2,4,model=False,portion='end') # port probs for kl divergence calculation 
            ses2_pps_model       = baseline_prob_dist(AQUA_visits2,4,model=True,portion='end') 
            
            across_ses_kl.append(cont_kl_divergence(lick_df2,lick_df3,ses2_pps,lim2,lim3,win_len,model=False))
            across_ses_kl_model.append(cont_kl_divergence(AQUA_visits2,AQUA_visits3,ses2_pps_model,lim2,lim3,win_len,model=True))

    if plot:
        plot_multises_kl(across_ses_kl,lim2,data='beh')
        plot_multises_kl(across_ses_kl_model,lim2,data='model')
        

def baseline_prob_dist(lick_df,ses_fraction,model=False,portion='end'):
    '''
    Compute a baseline port choice probability distribution from a fraction of a session.

    This distribution is used as the reference for KL divergence calculations. By default,
    it is computed from the final fraction of the session (e.g., the last quarter if
    `ses_fraction=4`), but can also be computed from the beginning of the session.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral or model-derived visit data. Must contain a 'port' column. If `model=False`,
        must also contain a 'unique_visit' column to filter first entries into visits.
    ses_fraction : int or float
        Fraction of the session to use (e.g., 4 uses the last or first quarter of the session).
    model : bool, optional
        If True, uses all rows in `lick_df['port']` (model-generated visits).
        If False, filters to unique visits using `unique_visit == 1`. Default is False.
    portion : {'end','beginning'}, optional
        Whether to use the final or initial fraction of the session. Default is 'end'.

    Returns
    -------
    port_probs : pandas.Series
        Probability of choosing each port (ports 1–6). Missing ports are assigned a small
        nonzero probability to avoid issues in downstream KL divergence calculations.
    '''
    if model:
        visits = lick_df['port']
    else:
        visits = lick_df[lick_df['unique_visit']==1]['port']
        
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


def cont_kl_divergence(lick_df2,lick_df3,port_probs,lim2,lim3,win_len=50,model=False):
    '''
    Calculate KL divergence from a reference port probability distribution.

    Computes choice distributions in sliding windows across the end of session 2
    and beginning of session 3, then calculates KL divergence between each windowed
    distribution and `port_probs`.

    Parameters
    ----------
    lick_df2 : pandas.DataFrame
        Session 2 behavioral or model visit data.
    lick_df3 : pandas.DataFrame
        Session 3 behavioral or model visit data.
    port_probs : pandas.Series
        Reference port probability distribution, typically computed from the stable
        portion of session 2.
    lim2 : int
        Number of visits from the end of session 2 to include.
    lim3 : int
        Number of visits from the beginning of session 3 to include.
    win_len : int, optional
        Sliding window length, in visits, used to compute local choice distributions.
        Default is 50.
    model : bool, optional
        If True, uses all rows in the input dataframes as model visits. If False,
        filters behavioral data to unique visits using `unique_visit == 1`.

    Returns
    -------
    kl_from_ses2 : numpy.ndarray
        KL divergence values for each sliding window relative to `port_probs`.
    '''
    if model: 
        visits2 = lick_df2['port'].iloc[-lim2:]
        visits3 = lick_df3['port'].iloc[:lim3]
    else:
        visits2     = lick_df2[lick_df2['unique_visit']==1]['port'].iloc[-lim2:]
        visits3     = lick_df3[lick_df3['unique_visit']==1]['port'].iloc[:lim3]
    visitsp     = pd.concat([visits2,visits3]).reset_index(drop=True)

    port_props_mov = [conv_to_prop(visitsp.iloc[i:win_len+i]) for i in range(len(visitsp))]
    mov_props      = pd.concat(port_props_mov,axis=1).transpose().dropna()
    kl_from_ses2   = np.array([entropy(mov_props.iloc[i,:],port_probs) for i in range(0,len(mov_props))])

    return kl_from_ses2 

def conv_to_prop(visits):
    '''
    Convert a sequence of port visits into a six-port probability distribution.
    
    Parameters
    ----------
    visits : array-like or pandas.Series
        Sequence of visited ports.
    
    Returns
    -------
    vis_props : pandas.Series
        Probability of visits to each port (ports 1–6). Ports not visited in the
        input are assigned probability 0.
    '''
    
    visits      = pd.Series(visits)
    missing_opt = np.array(list(set(np.array([1,2,3,4,5,6]))-set(visits.unique())))
    vis_props   = visits.value_counts().sort_index()/np.sum(visits.value_counts())
    if len(missing_opt)>0:
        for miss in missing_opt:
            vis_props = pd.concat([vis_props,pd.DataFrame([0],index=[miss])]).sort_index()
    return vis_props


### PLOTTING


def plot_multises_kl(across_ses_kl,lim2,data):
    '''
    Plot KL divergence around the session transition.

    Individual gray lines show mice and the thick black line shows the across-mouse
    mean. KL traces are smoothed separately before and after the session boundary
    so smoothing does not cross the transition. The x-axis is shifted so that the
    plotted range spans -200 to 300 visits, with the session transition marked by
    a dashed vertical line.

    Parameters
    ----------
    across_ses_kl : list of array-like
        KL divergence traces for each mouse.
    lim2 : int
        Boundary index separating the session 2 and session 3 portions of each KL trace.
    data : str
        Label indicating the data type being plotted, e.g. 'beh' or 'model'.

    Returns
    -------
    None
        Generates a matplotlib figure.
    '''
    smooth_kl = []
    for i,mouse in enumerate(across_ses_kl):
        smooth1 = pd.DataFrame(across_ses_kl[i][:lim2]).rolling(window=50,min_periods=2).mean()
        smooth2 = pd.DataFrame(across_ses_kl[i][lim2:]).rolling(window=50,min_periods=2).mean()
        smooth_kl.append(pd.concat([smooth1,smooth2]).reset_index(drop=True))
    smooth_kl = pd.concat(smooth_kl,axis=1)
        
    mean_kl = pd.DataFrame(across_ses_kl).transpose().mean(axis=1)
    mean_kl = pd.concat([mean_kl[:lim2].rolling(window=20,min_periods=2).mean(),mean_kl[lim2:].rolling(window=20,min_periods=2).mean()]).reset_index(drop=True)

    hfont   = {'fontname':'Arial'} 
    fig, ax = plt.subplots(figsize=(3,3))
    
    plot_len = len(mean_kl)
    x_vals   = np.arange(plot_len)-200

    ax.plot(x_vals,smooth_kl,color='gainsboro',linewidth=1.5)
    ax.plot(x_vals,mean_kl,color='black',linewidth=1.5)

    ax.axvline(x=lim2-200,color='darkgreen',linewidth=2,linestyle='--')

    ax.set_ylim([0,0.6])
    ax.set_yticks([0,1])
    ax.set_yticklabels([0,1],**hfont,fontsize=18)

    ax.set_xlim([-200,300])
    ax.set_xticks([-200,0,300])
    ax.set_xticklabels([-200,0,300],**hfont,fontsize=18)
    ax.set_xlabel('visit number',**hfont,fontsize=18)
    ax.set_ylabel('KL divergence (bits)',**hfont,fontsize=18)

    ax.spines[['right','top']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.spines['left'].set_position(('outward',10))
    ax.spines['bottom'].set_position(('outward',10))
    ax.tick_params(axis='both',which='both',labelsize=18,direction='in')
    
    fig.tight_layout()











