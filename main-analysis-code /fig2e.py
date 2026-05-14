import numpy as np
import seaborn as sns
import support_funcs as sf
import pandas as pd
import matplotlib.pyplot as plt
from itertools import permutations
from matplotlib.pyplot import cm

def MULTIivis(data_dict,ses_n=1,spatial_config=1,plot=True):
    """
    Compute transition metrics across mice for predefined port pairs.
    
    This function calculates:
    (1) `ivis` values for a set of port transition pairs, and
    (2) deviation from a null transition distribution (`target_trans_probs`).

    Mice are filtered by spatial configuration based on their interval structure.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary of loaded data.
    ses_n : int
        Session number to analyze.
    pair_list : list
        List of port-pair groupings (overwritten internally with predefined pairs).
    spatial_config : {1, 2}, default=1
        Select which spatial configuration of mice to include.
    plot : bool, default=True
        If True, plot transition probabilities relative to null.

    Returns
    -------
    pair_ivis : list
        List (per mouse) of `ivis` values for each pair set.
    """
    # all pairwise comparisons - hard coded 
    pair_list = [[(5,6),(1,2)],[(4,5),(1,3)],[(4,6),(2,3)],[(2,6)],[(2,4),(6,3)],[(2,5),(1,6)],[(1,4),(5,3)],[(4,3),(5,1)]]
    
    # extracting just interval task, non opto mice 
    subset_dict = sf.subset_mice(data_dict, task='conc', include_opto=False, config=1, region=None)
    
    pair_ivis, prob_from_null = [], []
    for mouse, mdata in subset_dict.items():
        is_config2 = mdata['conc']['b_meta']['intervals'][0][-1] == 1200
    
        if (spatial_config == 1 and is_config2) or (spatial_config == 2 and not is_config2):
            continue
    
        print(mouse)
    
        mlick_df = mdata['conc']['mlick_df']
        lick_df = mlick_df.loc[mlick_df['session_number'] == ses_n]
    
        if lick_df.empty:
            continue
    
        pair_ivis.append([ivis(pair, lick_df, split=False) for pair in pair_list])
        prob_from_null.append(target_trans_probs(lick_df, pair_list, seq_len=2))

    if plot:
        plot_prop_pairs(prob_from_null,spatial_config)
        
    return pair_ivis


def prop_ivis(pair_ivis):
    """
    Compute relative transition frequencies for each pair per mouse.

    For each mouse, this function converts counts of transitions (ivis) into
    proportions, so that the contributions from all pair sets sum to 1.

    Parameters
    ----------
    pair_ivis : list
        List (per mouse) of transition arrays for each pair set.

    Returns
    -------
    prop_pairs : list
        List (per mouse) of proportions for each pair set.
    """
    prop_pairs = []
    for mouse_ivis in pair_ivis:
        tot_ivis   = len(np.concatenate(mouse_ivis))
        prop_pairs.append([len(mouse_ivis[pair])/tot_ivis for pair in range(0,len(mouse_ivis))])
        
    return prop_pairs

def ivis(port_pair,lick_df,split=False):
    """
    Compute inter-visit intervals for specified port transitions.

    For a given set of port pairs, this function extracts the time between
    consecutive visits where the animal transitions between those ports.
    Intervals can be returned combined or separated by direction.

    Parameters
    ----------
    port_pair : list of tuple
        Port transitions to analyze, e.g. [(5, 6), (4, 5)].
    lick_df : pandas.DataFrame
        Lick/visit dataframe for one session.
    split : bool, default=False
        If True, return intervals separately for each direction
        (port1→port2 and port2→port1). If False, return combined intervals.

    Returns
    -------
    ivis : numpy.ndarray or tuple of numpy.ndarray
        Inter-visit intervals (seconds). If `split=True`, returns (ivis1, ivis2);
        otherwise returns a single concatenated array.
    """
    visits = lick_df.loc[lick_df['unique_visit']==1][['port','event_time']].reset_index(drop=True)
    # ivis1 = inter-visit intervals in one direction. ivis2 = in the other direction
    ivis1,ivis2 = [],[]
    for pair in port_pair: 
        ivis1.append(np.array([visits['event_time'].iloc[i+1] - visits['event_time'].iloc[i] for i in range(0,len(visits)-1) 
                  if ((visits['port'].iloc[i] == pair[0]) and (visits['port'].iloc[i+1] == pair[1]))]))
        ivis2.append(np.array([visits['event_time'].iloc[i+1] - visits['event_time'].iloc[i] for i in range(0,len(visits)-1) 
                  if ((visits['port'].iloc[i] == pair[1]) and (visits['port'].iloc[i+1] == pair[0]))]))

    if len(ivis1)==2:
        ivis1 = np.concatenate(ivis1)
        ivis2 = np.concatenate(ivis2)
    else:
        ivis1 = ivis1[0]
        ivis2 = ivis2[0]
    
    if split:
        return ivis1, ivis2
    else:
        return np.concatenate([ivis1,ivis2])
    
    
def sequencing(lick_df,seq_len=2):
    """
    Compute observed and null probabilities of port visit sequences.

    For all possible sequences of length `seq_len` (permutations of ports 1–6),
    this function calculates:
    - the observed probability of each sequence occurring in the data, and
    - a null probability assuming independent port visits.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Lick/visit dataframe for one session.
    seq_len : int, default=2
        Length of sequences to evaluate.

    Returns
    -------
    seq_df : pandas.DataFrame
        DataFrame with columns:
            - 'target' : sequence (list of ports)
            - 'prob'   : observed probability
            - 'null'   : expected probability under independence
            - 'diff'   : prob - null
    """
    
    port_visits = np.array(lick_df.loc[lick_df['unique_visit']==1]['port'])
    port_visits = [int(port) for port in port_visits]

    # generate target sequences through permutation
    comb = permutations([1,2,3,4,5,6], seq_len)
    
    all_seq_dict = []
    for target_sequence in comb:
        target_sequence = list(target_sequence)
        match = 0
        for i in range(0,len(port_visits)-1):
            if target_sequence == port_visits[i:i+len(target_sequence)]:
                match+=1
        seq_prob = match/(len(port_visits)-1)
    
        # calculate null (assumed independence): p1 x p2 x p3
        null_indep = np.prod([port_visits.count(port)/len(port_visits) for port in target_sequence])
            
        seq_dict = {'target': target_sequence,
                    'prob'  : seq_prob,
                    'null'  : null_indep}
        all_seq_dict.append(seq_dict)
        
    seq_df = pd.DataFrame(all_seq_dict)
    seq_df['diff'] = seq_df['prob']-seq_df['null']

    return(seq_df)


def target_trans_probs(lick_df, pair_list, seq_len=2):
    """
    Aggregate deviation-from-null probabilities for specified transition pairs.

    Uses `sequencing()` to compute observed and null probabilities for all
    sequences of length `seq_len`, then sums (observed − null) values for
    user-defined transition pairs. For each pair, both directions (A→B and B→A)
    are included.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Lick/visit dataframe for one session.
    pair_list : list
        List of pair sets. Each element is a list of port pairs (tuples),
        e.g. [[(5,6), (1,2)], [(4,6)]]. Each set is summed independently.
    seq_len : int, default=2
        Length of sequences used in probability calculation.

    Returns
    -------
    tot_diff : list
        List of summed (prob − null) values for each pair set.
    """
    seq_df = sequencing(lick_df,seq_len)

    tot_diff = []
    for pair_set in pair_list:
        diff_from_null = []

        for pair in pair_set:
            rev_pair = pair[::-1]

            val1 = seq_df.loc[seq_df['target'].isin([list(pair)]), 'diff']
            val2 = seq_df.loc[seq_df['target'].isin([list(rev_pair)]), 'diff']

            v1 = val1.iloc[0] if not val1.empty else 0
            v2 = val2.iloc[0] if not val2.empty else 0

            diff_from_null.append(v1 + v2)

        tot_diff.append(np.sum(diff_from_null))

    return tot_diff



### PLOTTING ###
    

def plot_prop_pairs(prop_pairs,spatial_config):
    """
    Plot transition probabilities relative to null across pair distances.

    Each pair set is assigned a fixed pairwise distance, and values are plotted
    across mice using violin plots with overlaid boxplots. This is intended to
    show whether certain classes of transitions occur more or less often than
    expected from the null distribution.

    Parameters
    ----------
    prop_pairs : list
        Per-mouse values for each predefined pair set.
    spatial_config : {1, 2}
        Spatial configuration to plot. If set to 2, a subset of mice is used.

    Returns
    -------
    None
    """
    distances = [14,18,22.8,42,56,65.5,70,72.2] # hard-coded, in order a --> h
    
    if spatial_config==2:
        prop_pairs = prop_pairs[-4:-1]

    # averaging 
    reorg_prop = [[item[i] for item in prop_pairs] for i in range(0,len(prop_pairs[0]))]
    
    # plotting 
    fig,ax = plt.subplots(figsize=(5.5,3.5))
    hfont  = {'fontname':'Arial'}
    vcolor = cm.viridis(np.linspace(0.9,0.1,num=len(prop_pairs[0])))
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    
    for i in range(0,8):
        sns.violinplot(y=reorg_prop[i],x=distances[i],ax=ax,color=vcolor[i],fill=False,width=2,split=True,inner=None,native_scale=True)
        sns.boxplot(y=reorg_prop[i],x=distances[i]+1,ax=ax,fill=False,width=2,showfliers=False,showcaps=False,color=vcolor[i],native_scale=True)
        
    ax.set_ylabel('p(transition) from null',**hfont,fontsize=18)
    ax.set_xlabel('pairwise distance',**hfont,fontsize=18)

    ax.set_ylim([-0.2,0.2])
    ax.set_xlim([12,distances[-1]+3])
    ax.set_yticks([-0.2,0,0.2])
    ax.set_yticklabels([-0.2,0,0.2])
    ax.set_xticks(distances)
    ax.set_xticklabels(['a','b','c','d','e','f','g','h'])

    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')

    ax.spines[['left','bottom']].set_position(('outward', 10))

    fig.tight_layout()
