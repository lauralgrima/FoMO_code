import math
import support_funcs as sf
import numpy as np
import pandas as pd
import travel as trav
import rescorla_wagner_sim as rw
import import_mat as im

def gen_observed(lick_df,bmeta):
    """
    Generate observed variables for regression from behavioral data.
    
    Filters to unique visits and returns a DataFrame containing choice, reward,
    photometry index, timing, session identity, and port rank (relative option
    value within each session).
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing at least 'unique_visit', 'port', 'rewarded',
        'photo_i', 'event_time', and 'session_number'.
    bmeta : object
        Behavioral metadata used to map ports to ranked options per session.
    
    Returns
    -------
    pandas.DataFrame
        Observations for each unique visit with columns:
        ['port', 'rewarded', 'photo_i', 'event_time',
         'session_number', 'port_ranks'].
    """
    observations = lick_df.loc[lick_df['unique_visit']==1][['port','rewarded','photo_i','event_time','session_number']].reset_index(drop=True)
    sessions = lick_df['session_number'].unique()
    ranks    = []
    for i in sessions:
        port_dict   = sf.port_2_rank(bmeta,i)
        port_visits = observations[observations['session_number']==i]['port']
        port_ranks  = [int(port_dict[vis]) for vis in port_visits]
        ranks.append(np.array(port_ranks))
    observations['port_ranks'] = np.concatenate(ranks)
    
    return observations


def gen_observed_AQUA(AQUA_beh_mouse,lick_df,bmeta,ses_n):
    """
    Generate observed variables for regression from AQUA behavior.
    
    Uses real behavioral timing/session structure, but replaces port choices,
    rewards, and port ranks with AQUA-derived values for the specified session.
    
    Parameters
    ----------
    AQUA_beh_mouse : dict
        AQUA behavioral data for a single mouse.
    lick_df : pandas.DataFrame
        Real behavioral data used for timing and photometry alignment.
    bmeta : object
        Behavioral metadata used for port-rank mapping.
    ses_n : list
        Session number list; the first session is used.
    
    Returns
    -------
    pandas.DataFrame
        Observations containing photo index, event time, session number, AQUA port,
        AQUA reward status, and AQUA port rank.
    """
    ldf_obs      = gen_observed(lick_df,bmeta)
    visits       = np.array(im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses_n[0],bmeta)['port'])
    port_ranks   = np.array(im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses_n[0],bmeta)['port_rank'])
    rews         = np.array(im.get_rews_AQUA(AQUA_beh_mouse['AQUA_rew_invis'],lick_df,ses_n[0])['rewarded'])

    observations               = ldf_obs[['photo_i','event_time','session_number']].copy()
    observations['port']       = visits
    observations['rewarded']   = rews
    observations['port_ranks'] = port_ranks
    
    return observations


def gen_predictors(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,data_dict,lick_df,vid_df,bmeta,ses_n,data_type='real',task='prob',import_travels=True,no_travel=True):
    """
    Generate behavioral predictor variables for regression analyses.
    
    Builds a predictor matrix for one session, including recent choice/reward
    history, interaction terms, matching-style reward history, Rescorla-Wagner
    variables, GLM-derived predictors, task timing variables, port rank, and
    movement-related predictors.
    
    Parameters
    ----------
    travels_filepath : str
        Path to saved travel-distance data.
    div_from_rand_filepath : str
        Path to divergence-from-random timing data.
    mat_GLM_path : str
        Path to GLM/AQUA model data.
    GLM_prob_filepath : str
        Path to probabilistic-task GLM data.
    data_dict : dict
        Mouse/session data used for importing GLM predictors.
    lick_df : pandas.DataFrame or list
        Behavioral data. For AQUA, use [lick_df, AQUA_beh_mouse].
    vid_df : pandas.DataFrame
        Video-derived data used for travel calculations.
    bmeta : dict-like
        Behavioral metadata for the animal/session.
    predictors : list
        Predictor names to return. Currently unused; the full predictor matrix is
        returned.
    ses_n : list
        Session number list; the first entry is used.
    data_type : {'real', 'sim', 'AQUA'}, default='real'
        Type of behavioral data.
    task : {'prob', 'conc'}, default='prob'
        Task type, used to choose GLM and divergence-from-random inputs.
    import_travels : bool, default=True
        Whether to load precomputed travel distances.
    no_travel : bool, default=True
        Whether to skip travel-distance calculation and fill with NaNs.
    
    Returns
    -------
    pandas.DataFrame
        Predictor matrix with behavioral, GLM-derived, timing, and movement
        predictors.
    """
    if data_type == 'AQUA':
        AQUA_beh_mouse = lick_df[1]
        lick_df        = lick_df[0]
    
    # dependent variables 
    if data_type == 'real':
        obs = gen_observed(lick_df,bmeta)
    elif data_type == 'sim':
        obs = lick_df 
    elif data_type == 'AQUA':
        obs = gen_observed_AQUA(AQUA_beh_mouse,lick_df,bmeta,ses_n)

    # whether n-x choice was at the same port
    n1_choice           = (obs['port'].shift(1) == obs['port']).astype(int).replace({0:-1,1:1})
    n1_choice           = n1_choice-n1_choice.mean()
    n1_choice.iloc[0]   = np.nan # change first one to nan as there is no previous
    n2_choice           = (obs['port'].shift(2) == obs['port']).astype(int).replace({0:-1,1:1})
    n2_choice           = n2_choice-n2_choice.mean()
    n2_choice.iloc[0:2] = np.nan
    n3_choice           = (obs['port'].shift(3) == obs['port']).astype(int).replace({0:-1,1:1})
    n3_choice           = n3_choice-n3_choice.mean()
    n3_choice.iloc[0:3] = np.nan
    
    # whether n-x choice was rewarded (same regardless of port). If wanting to change coding of variable, add .replace(0,-1) to end to change 0s to -1
    n1_rew = obs['rewarded'].shift(1).rename('n1_rew').replace({0:-1,1:1})
    n1_rew = n1_rew-n1_rew.mean()
    n2_rew = obs['rewarded'].shift(2).rename('n2_rew').replace({0:-1,1:1})
    n2_rew = n2_rew-n2_rew.mean()
    n3_rew = obs['rewarded'].shift(3).rename('n3_rew').replace({0:-1,1:1})
    n3_rew = n3_rew-n3_rew.mean()
    
    # interaction terms - for WSLS
    n1_int = (n1_rew*n1_choice)-(n1_rew*n1_choice).mean()
    n2_int = (n2_rew*n2_choice)-(n2_rew*n2_choice).mean()
    n3_int = (n3_rew*n3_choice)-(n3_rew*n3_choice).mean()
    
    # matching-type variables
    prop_rew_lo  = obs.groupby('port').cumsum()['rewarded']/obs.groupby('port').cumcount() # proportion of checks at that port that have been rewarded
    prop_rew_glo = obs.groupby('port').cumsum()['rewarded']/obs['rewarded'].cumsum() # what proportion of rewarded checks (of all rewarded checks across ports) have been at that port? 
    prop_rew_lo  = prop_rew_lo.replace(np.inf,np.nan)
    prop_rew_glo = prop_rew_glo.replace(np.inf,np.nan)
    
    # matching
    ratio_match = obs.groupby('port').cumsum()['rewarded']/(obs['rewarded'].cumsum()-obs.groupby('port').cumsum()['rewarded']) # ratio of rewards garnered at that port, i.e. 1/2-6 etc. 
    ratio_match = pd.DataFrame(sf.NormalizeNeg(ratio_match.replace(np.inf,np.nan)))
    ratio_match = ratio_match-ratio_match.mean()

    # cumulative reward by port, not normalised by number of visits (but normalised overall)
    cumu_rew_port = pd.DataFrame(sf.NormalizeNeg(obs.groupby('port').cumsum()['rewarded']))
    cumu_rew_port = cumu_rew_port-cumu_rew_port.mean()
    
    # cumulative reward overall
    cumu_rew = pd.DataFrame(sf.NormalizeNeg(obs['rewarded'].cumsum()))
    cumu_rew = cumu_rew-cumu_rew.mean()
    
    # rescorla wagner values 
    RWs   = rw.rescorla_wagner(lick_df,0.1,plot=False)
    value = np.array(sf.NormalizeNeg(np.array(RWs['vals'])))
    value = value-value.mean()
    pes   = np.array(sf.NormalizeNeg(np.array(RWs['pes'])))
    pes   = pes-pes.mean()
    
    # alpha, etc. - only for rewarded trials 
    if task == 'conc':
        aGLM_dfs,_   = im.import_GLMmat_data(mat_GLM_path,data_dict) # for all animals
    elif task == 'prob':
        aGLM_dfs = im.import_GLMmat_prob(GLM_prob_filepath,data_dict)
    mouse_GLM_df = [aGLM_dfs[i]  for i in range(len(aGLM_dfs)) if aGLM_dfs[i]['mouse'].iloc[0] == bmeta['animal_no']][0]

    # GLM values are for sessions 1 + 2 - limit 
    obs_nan = obs.dropna().reset_index(drop=True)
    if data_type == 'real':
        if ses_n[0] == 1:
            mouse_GLM_df_lim = mouse_GLM_df[:len(obs_nan[obs_nan['rewarded']==1])].reset_index(drop=True)
        elif ses_n[0] == 2:
            mouse_GLM_df_lim = mouse_GLM_df[len(mouse_GLM_df)-len(obs_nan[obs_nan['rewarded']==1]):].reset_index(drop=True)
        else:
            print('no')
    
        if not np.all(np.array(mouse_GLM_df_lim['port_id'])==np.array(obs_nan[obs_nan['rewarded']==1]['port'].astype(int))):
            print('mistmatch between GLM predictors + data')
        
        # rescale alpha, RPE variables - standardise 
        mouse_GLM_df_lim['alpha']    = sf.NormalizeNeg(mouse_GLM_df_lim['alpha'])-np.nanmean(sf.NormalizeNeg(mouse_GLM_df_lim['alpha']))
        mouse_GLM_df_lim['AQUA_rpe'] = sf.NormalizeNeg(mouse_GLM_df_lim['AQUA_rpe'])-np.nanmean(sf.NormalizeNeg(mouse_GLM_df_lim['AQUA_rpe']))
        if task == 'conc':
            mouse_GLM_df_lim['Q_rpe']    = sf.NormalizeNeg(mouse_GLM_df_lim['Q_rpe'])-np.nanmean(sf.NormalizeNeg(mouse_GLM_df_lim['Q_rpe']))
        
        # adding nans for unrewarded trials
        obs_nan['alpha'] = np.nan
        obs_nan.loc[obs_nan['rewarded']==1,'alpha'] = np.array(mouse_GLM_df_lim['alpha'])
        obs_nan['AQUA_rpe'] = np.nan
        obs_nan.loc[obs_nan['rewarded']==1,'AQUA_rpe'] = np.array(mouse_GLM_df_lim['AQUA_rpe'])
        if task == 'conc':
            obs_nan['Q_rpe'] = np.nan
            obs_nan.loc[obs_nan['rewarded']==1,'Q_rpe'] = np.array(mouse_GLM_df_lim['Q_rpe'])
        else:
            obs_nan['Q_rpe'] = np.nan
        
    else:
        obs_nan['alpha'] = np.nan
        obs_nan['Q_rpe'] = np.nan
        obs_nan['AQUA_rpe'] = np.nan
        
    # after random - 0 is before random, 1 is after 
    if task == 'conc':
        diverge_lats = sf.import_div_time(filepath=div_from_rand_filepath)
        try:
            mouse_lat    = diverge_lats[diverge_lats['Mouse']==bmeta['animal_no']]['Median latency (minutes)'].iloc[0] # checking if latency exists - horrible code
        except IndexError:
            mouse_lat    = np.nan
    elif task == 'prob':
        diverge_lats = {'dProb1' : 26.367225,
                        'dProb3' : 140.448767,
                        'gProb1' : 46.147242,
                        'gProb2' : 5.853892}
        mouse_lat = diverge_lats[bmeta['animal_no']]
        
    if not math.isnan(mouse_lat):
        after_random = (len(obs[obs['event_time']<(mouse_lat*60)])*[-0.5])+(len(obs[obs['event_time']>(mouse_lat*60)])*[0.5])
        after_random = sf.NormalizeNeg(np.array(after_random))
        after_random = after_random-np.nanmean(after_random)
    else:
        after_random = len(obs)*[-0.5] # i.e. never gets past random
        after_random = sf.NormalizeNeg(np.array(after_random))
        after_random = after_random-np.nanmean(after_random)

    # time since last check at that port
    ctime = obs.groupby('port')['event_time'].diff()
    ctime = sf.NormalizeNeg(ctime)
    ctime = ctime-np.nanmean(ctime)
    
    # time since last rewarded check at that port
    obs['ctime']   = obs.groupby('port')['event_time'].diff()
    time_since_rew = []
    for port in range(1,7):
        port_vis = obs.loc[obs['port']==port]
        port_vis = port_vis.copy()
        port_vis['group'] = port_vis['rewarded'].cumsum()
        port_vis['cumu_time'] = port_vis.groupby('group')['ctime'].cumsum().fillna(0)
        time_since_rew.append(port_vis['cumu_time'])
    rew_ctime = pd.Series(np.array(pd.DataFrame(time_since_rew).transpose().stack().sort_index()))
    rew_ctime = pd.Series(sf.NormalizeNeg(rew_ctime))
    rew_ctime = rew_ctime-rew_ctime.mean()

    # distance from previous port to current port (physical, not travelled)  
    distances = [[0,14,18,70,72.2,65.5],[14,0,22.8,56,65.5,42],[18,22.8,0,72.2,70,56],
                 [70,56,72.2,0,18,22.8],[72.2,65.5,70,18,0,14],[65.5,42,56,22.8,14,0]]
    dist = pd.Series([distances[int(obs['port'].iloc[i])-1][int(obs['port'].iloc[i+1])-1] for i in range(0,len(obs)-1)]).shift(1)
    dist = pd.Series(sf.NormalizeNeg(dist))
    dist = dist-dist.mean()
    
    # real distance travelled from previous port to current port 
    if data_type == 'real': # doesn't work for sim data 
        if not no_travel: #don't do travel at all
            if import_travels:
                saved_travels = import_saved_travels(csv_folder=travels_filepath,ses_n=[1,2])
                nvis          = len(lick_df[lick_df['unique_visit']==1])
                travel        = saved_travels[ses_n[0]-1][bmeta['animal_no']][:nvis].reset_index(drop=True)
                travel        = sf.NormalizeNeg(travel)
                travel        = travel-travel.mean()
            else: 
                travel = trav.dist_trav_between_visits(lick_df,bmeta,vid_df)
                travel = sf.NormalizeNeg(travel['path_length']).reset_index(drop=True)
        elif no_travel:
            travel = [np.nan]*(len(dist)+1)
    elif data_type == 'sim':
        travel = [np.nan]*(len(dist)+1)
    elif data_type == 'AQUA':
        travel = [np.nan]*(len(dist)+1)
        
    # port quality - flipping so that larger numbers are better, not worse (for DA regression interpretation)
    map_dict      = {1: 6, 2: 5, 3: 4, 4: 3, 5: 2, 6: 1}
    flipped_ranks = sf.NormalizeNeg(obs['port_ranks'].map(map_dict))-np.mean(sf.NormalizeNeg(obs['port_ranks'].map(map_dict)))
    
    # previous port
    prev_port = pd.DataFrame(flipped_ranks).shift(1)

    # create dataframe of predictors
    pred_df = pd.DataFrame(data={'rewarded'       :    obs['rewarded'].astype(int).replace({0:-0.5,1:0.5}),
                                 'port_ID'        :    obs['port'],
                                 'port_quality'   :    flipped_ranks,
                                 'prev_port_qual' :    prev_port.iloc[:,0],
                                 'n1_rew'         :    n1_rew, # same as stay/leave 
                                 'n2_rew'         :    n2_rew,
                                 'n3_rew'         :    n3_rew,
                                 'n1_choice'      :    n1_choice,
                                 'n2_choice'      :    n2_choice,
                                 'n3_choice'      :    n3_choice,
                                 'n1_int'         :    n1_int,
                                 'n2_int'         :    n2_int,
                                 'n3_int'         :    n3_int,
                                 'prop_rew_lo'    :    prop_rew_lo,
                                 'prop_rew_glo'   :    prop_rew_glo,
                                 'ratio_match'    :    ratio_match.iloc[:,0],
                                 'cumu_rew_port'  :    cumu_rew_port.iloc[:,0],
                                 'cumu_rew'       :    cumu_rew.iloc[:,0],
                                 'value'          :    value,
                                 'RW_pes'         :    pes,
                                 'alpha'          :    obs_nan['alpha'],
                                 'AQUA_rpe'       :    obs_nan['AQUA_rpe'],
                                 'Q_rpe'          :    obs_nan['Q_rpe'],
                                 'after_rand'     :    after_random,
                                 'time'           :    ctime,
                                 'time_s_rew'     :    rew_ctime,
                                 'dist'           :    dist,
                                 'travel'         :    travel})

    return pred_df


def gen_trav_predictor(data_dict,ses_n):
    """
    Generate travel-distance predictors across mice for selected sessions.
    
    Computes distance traveled between visits for each mouse/session and returns
    the normalized travel predictor in a wide DataFrame, with one column per mouse.
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral and video data.
    ses_n : list
        Sessions to include, typically a single session such as [1].
    
    Returns
    -------
    pandas.DataFrame
        Normalized travel-distance predictors, with mice as columns.
    """
    travels,mouseys = [],[]
    for i, mouse in enumerate(list(data_dict.keys())):
        print(mouse)
        
        tot_ses   = data_dict[mouse]['b_meta']['nsessions'] # total number of sessions for given mouse 
        ses_count = list(range(1,tot_ses+1))
        sessions  = list(set(ses_count).intersection(ses_n))
        if sessions: # if not an empty list - ie if the mouse has done sessions in the range relevant, as given by ses_n
            mouseys.append(mouse)
            for ses in sessions: 
                lick_df,_,vid_df,bmeta,_ = sf.extract_data(data_dict,mouse,ses) 
                travel                   = trav.dist_trav_between_visits(lick_df,bmeta,vid_df)
                travel                   = sf.NormalizeData(travel['path_length']).reset_index(drop=True)
                travels.append(travel)

    # reformat 
    travels         = pd.concat(travels,axis=1)
    travels.columns = mouseys

    return travels 
    
def import_saved_travels(csv_folder,ses_n=[1,2]):
    """
    Load precomputed travel-distance predictors from CSV files.
    
    Reads saved travel-distance data for specified sessions and returns them
    as a list of DataFrames, one per session.
    
    Parameters
    ----------
    csv_folder : str
        Directory containing saved travel CSV files. Files should be named
        'travels<session>.csv'.
    ses_n : list, default=[1, 2]
        Sessions to load.
    
    Returns
    -------
    list of pandas.DataFrame
        Travel-distance data for each requested session.
    """
    saved_travels = []
    for ses in ses_n:
        saved_travels.append(pd.read_csv(csv_folder+'travels'+str(ses)+'.csv',index_col=0))
    
    return saved_travels


