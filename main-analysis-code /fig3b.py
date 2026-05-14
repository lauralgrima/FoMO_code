import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import import_mat as im
import support_funcs as sf


def eg_transition_matrix(data_dict,mat_GLM_path,mouse='6PG6',port_or_rank='rank',ses_n=1,plot=True):
    """
    Calculate example transition matrices for one mouse and one session.

    Parameters:
        data_dict: dict
            Full behavioral data dictionary.

        mat_GLM_path: str
            Path to AQUA/GLM .mat data.

        mouse: str
            Mouse ID.

        port_or_rank: str
            'port' to calculate transitions between physical ports.
            'rank' to calculate transitions between port ranks.

        ses_n: int
            Session number.

        plot: bool
            If True, plot real and AQUA transition matrices.

    Returns:
        mouse_trans_mat: pd.DataFrame
            Transition matrix from real behavior.

        AQUA_trans_mat: pd.DataFrame
            Transition matrix from AQUA/model behavior.
    """
    
    lick_df  = data_dict[mouse]['conc']['mlick_df']
    bmeta    = data_dict[mouse]['conc']['b_meta']
    
    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    _,AQUA_beh_NAc  = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS  = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_beh        = {**AQUA_beh_DMS,**AQUA_beh_NAc}

    # calculate transition matrices - real and model behavior 
    mouse_trans_mat   = transition_matrix(lick_df,bmeta,ses_n,real_or_sim='real')
    
    AQUA_visits = im.get_visits_AQUA(AQUA_beh[mouse]['AQUA_vis_invis'],lick_df,ses_n,bmeta)['port']
    if port_or_rank == 'rank':
        AQUA_visits = AQUA_visits.map(sf.port_2_rank(bmeta,ses_n))
    AQUA_trans_mat = transition_matrix(AQUA_visits,bmeta,ses_n,real_or_sim='AQUA')
    
    if plot: 
        plot_ind_trans_mat(mouse_trans_mat)
        plot_ind_trans_mat(AQUA_trans_mat)


def MULTItransition_matrix(data_dict,mat_GLM_path,ses_n=1,port_or_rank='rank'):
    '''
    Calculate transition matrices for real and AQUA data across mice for one session.

    Notes:
        Currently for a single session only.
        Transition matrices sum to 1 across all 36 transitions, not per row.

    Parameters:
        data_dict: dict
            Full behavioral data dictionary.

        mat_GLM_path: str
            Path to AQUA/GLM .mat data.

            Session number to analyze.

        port_or_rank: str
            'port' to calculate transitions between physical ports.
            'rank' to calculate transitions between port ranks.

    Returns:
        all_trans_mat: list
            List of real-data transition matrices, one per mouse.

        all_AQUA_trans_mat: list
            List of AQUA transition matrices, one per mouse.
    '''
    # extracting just interval task, non opto mice, config 1
    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    subset_dict     = {**subset_dict_NAc,**subset_dict_DMS}

    # importing AQUA model data 
    _,AQUA_beh_NAc = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_beh = {**AQUA_beh_DMS,**AQUA_beh_NAc}
    
    all_trans_mat,all_AQUA_trans_mat,all_bmeta = [],[],[]
    for mouse in subset_dict:
        print(mouse)
        
        tot_ses = data_dict[mouse]['conc']['b_meta']['nsessions']
        if ses_n <= tot_ses:
            lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(subset_dict,mouse,ses_n)
            trans_mat   = transition_matrix(lick_df,bmeta,ses_n,real_or_sim='real')
            AQUA_visits = im.get_visits_AQUA(AQUA_beh[mouse]['AQUA_vis_invis'],lick_df,ses_n,bmeta)['port']
            if port_or_rank == 'rank':
                AQUA_visits = AQUA_visits.map(sf.port_2_rank(bmeta,ses_n))
            AQUA_trans_mat = transition_matrix(AQUA_visits,bmeta,ses_n,real_or_sim='AQUA')
        all_bmeta.append(bmeta)
        all_trans_mat.append(pd.DataFrame(trans_mat))
        all_AQUA_trans_mat.append(pd.DataFrame(AQUA_trans_mat))

    return all_trans_mat, all_AQUA_trans_mat


def transition_matrix(lick_df,bmeta,n_ses,port_or_rank='port',real_or_sim='real'):
    '''
    Calculate probability of transitioning from one port to another for a given session.

    Returns:
        trans_matrix: pandas DataFrame (6x6)
            Rows = current port (1–6)
            Columns = next port (1–6)
            Values = probability of transitioning

    Variables:
        port_or_rank: 
            'port' → transitions between physical ports
            'rank' → transitions between ranked ports
        real_or_sim: 
            'real' → real data
            'sim'  → simulated data (port_or_rank ignored)
            'AQUA' → AQUA data (lick_df already formatted as sequence)

    Notes:
        FOR SIM DATA:
            lick_df should already be a sequence of port visits
    '''
    
    if real_or_sim=='real':
        port_visits  = lick_df.loc[lick_df['unique_visit']==1]['port']
        port_visits  = [int(port_visits.iloc[i]) for i in range(0,len(port_visits))] 
        if port_or_rank == 'rank': # convert to rank 
            port_dict   = sf.port_2_rank(bmeta,n_ses)
            port_visits = [int(port_dict[vis]) for vis in port_visits]
        port_visits = pd.Series(port_visits) # convert back to pandas 
    elif real_or_sim == 'AQUA':
        port_visits = lick_df
        
    all_counts   = []
    for port_n in range(1,7):
        next_port      = [port_visits.iloc[i+1] for i,visit in enumerate(port_visits.iloc[:-1]) if visit == port_n]
        unique, counts = np.unique(next_port, return_counts=True)
        missing_i      = ((np.setdiff1d([1,2,3,4,5,6],unique))-1) # if a specific transition was not made 
        for i in missing_i:
            counts = np.insert(counts,i,0)
        all_counts.append(counts)
    # each row is by port - e.g. row 1 = transitions from port 1 to other ports (including itself)
    trans_matrix = all_counts/np.sum(np.concatenate(all_counts))

    return pd.DataFrame(trans_matrix,range(1,7),range(1,7))
        

def plot_ind_trans_mat(trans_mat):
    '''
    Plot individual transition matrix (6x6) as a heatmap.

    Parameters:
        trans_mat (pd.DataFrame or array-like): 
            6x6 transition probability matrix (rows = current, cols = next)
    '''
    fig, ax = plt.subplots(figsize=(3, 3))

    sns.heatmap(
        trans_mat,
        ax=ax,
        vmin=0,
        vmax=0.2,
        cmap='Greens',
        cbar=False,
        square=True,
        linewidths=0.5,
        linecolor='white'
    )

    # ticks and labels
    ax.xaxis.tick_top()
    ax.set_xlabel('Next port', labelpad=10)
    ax.set_ylabel('Current port')

    ax.set_xticklabels(range(1, 7))
    ax.set_yticklabels(range(1, 7), rotation=0)

    # cleaner look
    ax.tick_params(left=False, bottom=False, top=False)

    fig.tight_layout()