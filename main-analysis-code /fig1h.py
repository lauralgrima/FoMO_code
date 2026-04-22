import scipy
import pandas as pd 
import numpy as np
import seaborn as sns
import support_funcs as sf
import strategy_benchmarks as sb
import matplotlib.pyplot as plt



def MULTIkl_diverge_ratio(data_dict,ses_n=1,final_policy='IO',consec_counts=100,nrand_runs=100,plot_ind=False,plot=True):
    """
    Compute KL-ratio threshold metrics across mice for one session.

    Subsets to non-opto mice with the 'conc' task, runs `kl_diverge_ratio()`
    for each mouse, and optionally plots the across-mouse distribution.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary of loaded data.
    ses_n : int, default=1
        Session number to analyze.
    final_policy : {'IO', 'matching'}, default='IO'
        Reference policy used in the KL ratio.
    consec_counts : int, default=100
        Number of consecutive samples above threshold required.
    nrand_runs : int, default=100
        Number of random simulations per mouse.
    plot_ind : bool, default=False
        If True, plot per-mouse KL-ratio examples.
    plot : bool, default=True
        If True, plot summary results across mice.

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

    med_nrew_IOs = []
    for mouse in subset_dict:
        print(mouse)
        lick_df,_,_,bmeta,_ = sf.extract_data(subset_dict,mouse,ses_n)
        med_nrew_IOs.append(kl_diverge_ratio(lick_df,bmeta,ses_n,final_policy,nrand_runs,consec_counts,plot=plot_ind))
    
        if plot: # plot one session at a time
            plot_multi_OIO(np.array(list(filter(None,med_nrew_IOs))),ylabel='n. rewards to ratio > 1')


def kl_diverge_ratio(lick_df,bmeta,ses_n,final_policy='IO',nrand_runs=10,consec_counts=100,plot=True):
    """
    Estimate when behavior becomes more like a target policy than random.

    Computes visit distributions for the animal, random policy, and a target
    policy ('IO' or 'matching'), then calculates the ratio of KL divergences
    from animal-to-random over animal-to-target.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Lick/visit dataframe for one session.
    bmeta : dict
        Behavioral metadata used for simulation.
    ses_n : int
        Session number.
    final_policy : {'IO', 'matching'}, default='IO'
        Target policy used in the KL-ratio calculation.
    nrand_runs : int, default=10
        Number of random simulations to average over.
    consec_counts : int, default=100
        Number of consecutive samples with ratio > 1 required.
    plot : bool, default=True
        If True, plot a single-session KL-ratio example.

    Returns
    -------
    med_nrew_IO : float or list
        Median number of rewards at which the KL ratio crosses threshold.
        Returns an empty list if threshold is not reached reliably.
    """
    # mouse pvisits
    visits    = lick_df[lick_df['unique_visit']==1][['port','rewarded']]
    cumcounts = pd.concat([visits.groupby('port').cumcount().add(1).where(visits['port'].eq(port)).ffill().fillna(0) for port in range(1,7)],axis=1)
    cum_probs = pd.DataFrame([np.array([cumcounts.iloc[j,i]/np.sum(cumcounts.iloc[j,:]) for i in range(0,6)]) for j in range(len(cumcounts))])

    # random pvisits
    visit_times     = np.array(lick_df.loc[lick_df['unique_visit']==1]['event_time'])
    rand_probs_runs = []
    for i in range(nrand_runs):
        print(i)
        rand_sim              = sb.simulate_strategies(visit_times,bmeta,ses_n,stagger=False,strategy='random')[['port_list','outcomes','rewards_col']]
        rand_sim['port_list'] = rand_sim['port_list']+1
        rand_cumcounts        = pd.concat([rand_sim.groupby('port_list').cumcount().add(1).where(rand_sim['port_list'].eq(port)).ffill().fillna(0) for port in range(1,7)],axis=1)
        rand_probs            = pd.DataFrame([np.array([rand_cumcounts.iloc[j,i]/np.sum(rand_cumcounts.iloc[j,:]) for i in range(0,6)]) for j in range(0,len(rand_cumcounts))])
        rand_probs_runs.append(rand_probs)

    # ideal observer pvisits 
    IO_sim = sb.simulate_strategies(visit_times,bmeta,ses_n,stagger=True,strategy='io')[['port_list','outcomes','rewards_col']]
    IO_sim['port_list'] = np.array([IO_sim['port_list'].iloc[i][0] if isinstance(IO_sim['port_list'].iloc[i],list) else IO_sim['port_list'].iloc[i] for i in range(0,len(IO_sim['port_list']))])+1
    IO_cumcounts        = pd.concat([IO_sim.groupby('port_list').cumcount().add(1).where(IO_sim['port_list'].eq(port)).ffill().fillna(0) for port in range(1,7)],axis=1)
    IO_probs            = pd.DataFrame([np.array([IO_cumcounts.iloc[j,i]/np.sum(IO_cumcounts.iloc[j,:]) for i in range(0,6)]) for j in range(0,len(IO_cumcounts))])

    # matching pvisits
    cum_rew         = pd.concat([visits.groupby('port').cumsum().where(visits['port'].eq(port)).ffill().fillna(0) for port in range(1,7)],axis=1)
    cum_rew.columns = [1,2,3,4,5,6]

    epsilon        = 0.0000000001 # to avoid dividing by zero
    row_sums       = np.sum(cum_rew.values, axis=1, keepdims=True)  # Sum across columns (per row)
    denominators   = row_sums - (cum_rew.values + epsilon)     
    matching_probs = pd.DataFrame(abs(cum_rew.values / denominators)) # element-wise division

    # remove first row for matching
    cum_probs         = cum_probs[1:]
    matching_probs    = matching_probs[1:]
    IO_probs          = IO_probs[1:]
    
    # replace zero with very small value
    matching_probs = matching_probs.replace(0,0.00001)
    cum_probs      = cum_probs.replace(0,0.00001)
    IO_probs       = IO_probs.replace(0,0.00001)
    
    # kl divergence calculation
    kl_from_match = np.array([scipy.stats.entropy(cum_probs.iloc[i,:],matching_probs.iloc[i,:]) for i in range(0,len(cum_probs))])
    kl_from_IO    = np.array([scipy.stats.entropy(cum_probs.iloc[i,:],IO_probs.iloc[i,:]) for i in range(0,len(IO_probs))])
    
    if len(cum_probs)>len(rand_probs):
        cum_probs      = cum_probs[:-1]
        matching_probs = matching_probs[:-1]
        
    if len(cum_probs)>len(IO_probs):
        cum_probs = cum_probs[:-1]

    # dealing with random as it has multiple runs
    nrew_IOs,nvis_IOs = [],[]
    for arr in rand_probs_runs:
        if len(arr)>len(cum_probs):
            arr          = arr[1:]
        arr          = arr.replace(0,0.00001)
        kl_from_rand = np.array([scipy.stats.entropy(cum_probs.iloc[i,:],arr.iloc[i,:]) for i in range(0,len(cum_probs))])
        
        if final_policy == 'IO':
            kl_ratio = (kl_from_rand/kl_from_IO)[50:]
        elif final_policy == 'matching':
            kl_ratio = (kl_from_rand/kl_from_match)[50:]
        
        iover_thresh = thresh_ratio(kl_ratio,consec_counts)

        # use iover_thresh to calculate n visits and n rewards when animal is more like OIO
        # probably should be in separate function in future 
        if isinstance(iover_thresh,np.int64):
            nvis_IOs.append(iover_thresh+50) # right now always ignoring first 50
            nrew_IOs.append(visits.iloc[:iover_thresh+50]['rewarded'].sum())
        else:
            nvis_IOs.append([])
            nrew_IOs.append([])
            
    if not any(not sublist for sublist in nrew_IOs):
        med_nrew_IO = np.median(nrew_IOs)
    else:
        empty_count = sum(1 for item in nrew_IOs if isinstance(item, list) and len(item) == 0)
        if empty_count>8:
            print('rubbish behaviour, be wary')
            med_nrew_IO = []
        else:
            med_nrew_IO = np.median(nrew_IOs)

    if plot:
        # plotting a single example 
        plot_eg_kl_ratio(kl_ratio,iover_thresh,[],bmeta,smooth=True)
        
    return med_nrew_IO


def thresh_ratio(kl_ratio,consec_counts=100):
    """
    Find the first index where KL ratio stays above 1 long enough.

    Parameters
    ----------
    kl_ratio : array-like
        KL-divergence ratio over visits.
    consec_counts : int, default=100
        Number of consecutive samples above 1 required.

    Returns
    -------
    iover_thresh : int or list
        First threshold-crossing index, or an empty list if none is found.
    """
    diff_rat     = np.diff(np.where(kl_ratio>1)[0])
    if diff_rat.size > 0:
        diff_mask    = diff_rat == 1 
        conv         = np.convolve(diff_mask, np.ones(consec_counts, dtype=int), mode='valid')
        iover_thresh = np.sum(diff_rat[:np.argmax(conv >= consec_counts)])+np.where(kl_ratio>1)[0][0]
    else:
        iover_thresh = []
    
    return iover_thresh
       

def plot_eg_kl_ratio(kl_ratio,iover_thresh,inew_ses,bmeta,smooth=True):
    """
    Plot KL-divergence ratio for a single session.

    Parameters
    ----------
    kl_ratio : array-like
        KL-divergence ratio over visits.
    iover_thresh : int or list
        Threshold-crossing index to mark on the plot.
    inew_ses : int or list
        Optional session-change index to mark on the plot.
    bmeta : dict
        Behavioral metadata, used here for the plot title.
    smooth : bool, default=True
        If True, smooth the ratio with a moving average.

    Returns
    -------
    None
    """
    hfont = {'fontname':'Arial'} 
    fig, ax = plt.subplots(figsize=(3.5,3.5))
    
    if smooth:
        kl_ratio_smooth = np.convolve(kl_ratio, np.ones(10)/10, mode='valid')
        ax.plot(kl_ratio_smooth,color='black')
    else:
        ax.plot(kl_ratio,color='black')
        
    if iover_thresh:
        ax.axvline(iover_thresh,color='darkgrey',linestyle='--',linewidth=2)
        
    if inew_ses:
        ax.axvline(inew_ses,color='purple',linestyle='--',linewidth=2)
        
    ax.axhline(y=1,color='darkgrey',linewidth=2)
    
    ax.set_yscale('log')

  #  ax.set_ylim([0,10])
    ax.set_xlim([0,len(kl_ratio)])
    ax.set_xticks([0,len(kl_ratio)])
    ax.set_yticks([1,10])
    ax.set_xticklabels([0,len(kl_ratio)],**hfont,fontsize=18)
    ax.set_yticklabels([1,10],**hfont,fontsize=18)
    ax.set_xlabel('visit number', **hfont,fontsize=18)
    ax.set_ylabel('KL divergence ratio (log)', **hfont,fontsize=18)
    ax.spines[['right','top']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=14,direction='in')
    
    ax.set_title(bmeta['animal_no'])
    
    fig.tight_layout()
    plt.show()
    
    
def plot_multi_OIO(values,ylabel):
    """
    Plot a summary distribution for a single across-mouse metric.

    Parameters
    ----------
    values : array-like
        Values to plot.
    ylabel : str
        Y-axis label.
    title : str
        Plot title or label identifier.

    Returns
    -------
    None
    """
    
    # plotting time to diverge from random 
    fig,ax = plt.subplots(figsize=(3,3))
    hfont  = {'fontname':'Arial'} 
    
    dot_colors  = len(values)*['grey']
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    sns.boxplot(y=values,x=0,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
    sns.stripplot(y=values,x=0,edgecolor=dot_colors,facecolor='white',linewidth=2)
 
    ax.set_ylabel(ylabel,**hfont,fontsize=18)
    
    ax.set_yticks([0,50,100,150,200])
    ax.set_yticklabels([0,50,100,150,200],fontsize=18)
    ax.set_xlim([-1,5])
    ax.set_ylim([0,250])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()
    plt.show()