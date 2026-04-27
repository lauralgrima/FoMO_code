import numpy as np
import pandas as pd
import support_funcs as sf
#import analysis.travel_and_paths as trav
#import analysis.rescorla_wagner as rw
import import_mat as im
from sklearn.preprocessing import StandardScaler

# generating predictor matrix and observation matrix for regressions (behaviour and photometry)
# general approach: create large predictor df that has all predictors, then in individual functions for running regressions reduce down to trial types needed 

def gen_observed(lick_df,bmeta):
    '''
    Dependent variable for regression.
    Returns:
        observed_choices   :   port chosen, 1 - 6
        observed_rews      :   whether choice was rewarded (1 or 0)
        observed_photo     :   photometry indices corresponding to choices 
        observed_time      :   time of choice 
        session number     :   session number (for splitting by session)
        port_ranks         :   rank of chosen option, 1 - 6 
    '''
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
    '''
    Observations but for AQUA format. This is for behaviour - comparing mouse and AQUA in behaviour regression. 
    AQUA_beh_mouse should be the dictionary for a single animal. 
    '''
    ldf_obs      = gen_observed(lick_df,bmeta)
    visits       = np.array(im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses_n,bmeta)['port'])
    port_ranks   = np.array(im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses_n,bmeta)['port_rank'])
    rews         = np.array(im.get_rews_AQUA(AQUA_beh_mouse['AQUA_rew_invis'],lick_df,ses_n)['rewarded'])

    observations               = ldf_obs[['photo_i','event_time','session_number']].copy()
    observations['port']       = visits
    observations['rewarded']   = rews
    observations['port_ranks'] = port_ranks
    
    return observations


def gen_predictors(data_dict,lick_df,vid_df,bmeta,ses_n,predictors,data_type='real',scale=True):
    '''
    Generates a predictor matrix for all variables relevant to behaviour and photometry regressions. 
    Only for a single session (determined by ses_n). data_dict is required as an argument just for importing one variable. 
    Note that scaling happens at the end and is z-scoring (mean of 0, SD of 1)
    
    Args:
        data_type      :     set to 'sim' for sim data - travel will just be NaN, 'real' for mouse bhehaviour, 'AQUA' for AQUA simulation specifically. 
        import_travels :     set to True to import previously calculated travels from csv (speeds things up)
        scale          :     set to True to z-score all predictors 
    
    Predictors:
        n1_rew,n2_rew,n3_rew            :   whether n-x choice was rewarded (regardless of port)
        n1_choice,n2_choice,n3_choice   :   whether n-x choice was at the same port 
        prop_rew_lo, prop_rew_glo       :   probability of reward at option locally and globally (proportion of rewarded visits at an individual option, and proportion of rewarded checks that have been at that option)
        prew                            :   probability of reward at all options over some window (currently hard-coded at 500)
        

    Interaction terms: 
        n1_interaction, n2_interaction  :  reward x choice, up to 3 choices back
    
    '''
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
    n1_choice.iloc[0]   = np.nan # change first one to nan as there is no previous
    n2_choice           = (obs['port'].shift(2) == obs['port']).astype(int).replace({0:-1,1:1})
    n2_choice.iloc[0:2] = np.nan
    n3_choice           = (obs['port'].shift(3) == obs['port']).astype(int).replace({0:-1,1:1})
    n3_choice.iloc[0:3] = np.nan
    
    # whether n-x choice was rewarded (same regardless of port). If wanting to change coding of variable, add .replace(0,-1) to end to change 0s to -1
    n1_rew = obs['rewarded'].shift(1).rename('n1_rew').replace({0:-1,1:1})
    n2_rew = obs['rewarded'].shift(2).rename('n2_rew').replace({0:-1,1:1})
    n3_rew = obs['rewarded'].shift(3).rename('n3_rew').replace({0:-1,1:1})
    
    # interaction terms - for WSLS
    n1_int = (n1_rew*n1_choice)-(n1_rew*n1_choice).mean()
    n2_int = (n2_rew*n2_choice)-(n2_rew*n2_choice).mean()
    n3_int = (n3_rew*n3_choice)-(n3_rew*n3_choice).mean()
    
    # matching-type variables
    prop_rew_lo  = obs.groupby('port').cumsum()['rewarded']/obs.groupby('port').cumcount() # cumulative proportion of checks at that port that have been rewarded
    prop_rew_glo = obs.groupby('port').cumsum()['rewarded']/obs['rewarded'].cumsum() # what proportion of rewarded checks (of all rewarded checks across ports) have been at that port? cumulative 
    prop_rew_lo  = prop_rew_lo.replace(np.inf,np.nan)
    prop_rew_glo = prop_rew_glo.replace(np.inf,np.nan)
    
    # probability of reward 
    prew = pd.Series(obs['rewarded']).rolling(window=500,center=True,min_periods=1).mean()
    
    # matching
    ratio_match = obs.groupby('port').cumsum()['rewarded']/(obs['rewarded'].cumsum()-obs.groupby('port').cumsum()['rewarded']) # ratio of rewards garnered at that port, i.e. 1/2-6 etc. 
    ratio_match = ratio_match.replace([np.inf,-np.inf],np.nan)

    # cumulative reward by port, not normalised by number of visits (but normalised overall)
    cumu_rew_port = obs.groupby('port').cumsum()['rewarded']
    
    # cumulative reward overall
    cumu_rew = obs['rewarded'].cumsum()

    # time since last check at that port
    ctime = obs.groupby('port')['event_time'].diff()
    
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

    # distance from previous port to current port (physical, not travelled)  
    distances = [[0,14,18,70,72.2,65.5],[14,0,22.8,56,65.5,42],[18,22.8,0,72.2,70,56],
                 [70,56,72.2,0,18,22.8],[72.2,65.5,70,18,0,14],[65.5,42,56,22.8,14,0]]
    dist = pd.Series([distances[int(obs['port'].iloc[i])-1][int(obs['port'].iloc[i+1])-1] for i in range(0,len(obs)-1)]).shift(1)
    dist.at[len(dist)] = np.nan

    # previous port
    prev_port = obs['port_ranks'].shift(1)

    # create dataframe of predictors
    pred_df = pd.DataFrame(data={'rewarded'       :    obs['rewarded'].astype(int).replace({0:-0.5,1:0.5}),
                                 'port_ID'        :    obs['port'],
                                 'port_quality'   :    obs['port_ranks'],
                                 'prev_port_qual' :    prev_port,
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
                                 'prew'           :    prew,
                                 'ratio_match'    :    ratio_match,
                                 'cumu_rew_port'  :    cumu_rew_port,
                                 'cumu_rew'       :    cumu_rew,
                                 'time'           :    ctime,
                                 'time_s_rew'     :    rew_ctime,
                                 'dist'           :    dist})

    if scale:
        scaler      = StandardScaler()
        scaled_data = pd.DataFrame(scaler.fit_transform(pred_df))
        scaled_data.columns = pred_df.columns
        pred_df     = scaled_data 

    return pred_df[predictors]



