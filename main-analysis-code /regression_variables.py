import numpy as np
import pandas as pd
import support_funcs as sf
import analysis.travel_and_paths as trav
import analysis.rescorla_wagner as rw
import preprocessing.import_mat as im
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
    visits       = np.array(im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses_n[0],bmeta)['port'])
    port_ranks   = np.array(im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses_n[0],bmeta)['port_rank'])
    rews         = np.array(im.get_rews_AQUA(AQUA_beh_mouse['AQUA_rew_invis'],lick_df,ses_n[0])['rewarded'])

    observations               = ldf_obs[['photo_i','event_time','session_number']].copy()
    observations['port']       = visits
    observations['rewarded']   = rews
    observations['port_ranks'] = port_ranks
    
    return observations






