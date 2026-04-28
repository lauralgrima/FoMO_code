import numpy as np
import pandas as pd
import preprocessing.support_funcs as sf
from scipy.io import loadmat
from itertools import chain


### MAIN IMPORTING FUNCS

def import_GLMmat_data(path,data_dict,region='NAc'):
    '''
    regions can be either ['NAc','DMS'] or just ['NAc']
    Path should be folder which contains files (files are separate for DMS and NAc, but will be combined). Path doesn't have to have slash at end
    This function is specifically for MATLAB data to do with values for dopamine GLM
    Returns a list of dataframes, one df per animal. Each df contains:
        mouse    :    mouse ID
        alpha    :    np.array of alpha values, one value per reward 
        AQUA_rpe :    np.array of RPE values from AQUA model, one value per reward
        Q_rpe    :    np.array of RPE values from q learning model, one value per reward 
        port_id  :    np.array of actual port sampled, one value per rewarded visit 
    Region is not explicitly 'encoded' in any of these variables (but should be fine if selecting by mouse ID)
    '''
    
    if region == 'NAc':
        full_path = path+'GLM_export_v7_final_NAc.mat'

    elif region == 'DMS':
        full_path = path+'GLM_export_v7_final_DMS.mat'

    mat_data = loadmat(full_path)
    
    # reformat
    data = mat_data['GLM_export'][0]
    
    GLM_dfs = []
    AQUA_dicts = {}
    for i in range(len(data)):
        # variables for GLM analyses 
        GLM_dict               = {}
        GLM_dict['mouse']      = data[i]['anim_id'][0][:6]
        GLM_dict['alpha']      = np.concatenate(data[i]['alpha'])
        GLM_dict['AQUA_rpe']   = np.concatenate(data[i]['AQUA_rpe'])
        GLM_dict['Q_rpe']      = np.concatenate(data[i]['Q_rpe'])
        GLM_dict['port_id']    = np.concatenate(data[i]['port_id'])
        GLM_df                 = pd.DataFrame(GLM_dict)
        GLM_dfs.append(GLM_df)
        
        # behavioural output of AQUA model 
        AQUA_dict                   = {}
        AQUA_dict['mouse']          = data[i]['anim_id'][0][:6]
        AQUA_dict['AQUA_rew']       = pd.DataFrame(data[i]['aqua_model']['rewards_for_LL'][0][0]).transpose()
        AQUA_dict['AQUA_vis']       = pd.DataFrame(data[i]['aqua_model']['visits_for_LL'][0][0]).transpose()
        AQUA_dict['AQUA_rew_invis'] = pd.DataFrame(data[i]['aqua_model']['rewards_invis'][0][0]).transpose()
        AQUA_dict['AQUA_vis_invis'] = pd.DataFrame(data[i]['aqua_model']['visits_invis'][0][0]).transpose()
        AQUA_dict['r2']             = np.concatenate(data[i]['best_r2'])
        # adding useful stuff from data_dict 
        try:
            AQUA_dict['b_meta']    = data_dict[AQUA_dict['mouse']]['conc']['b_meta']
            ecounts     = int(np.ceil(len(AQUA_dict['AQUA_vis'])/AQUA_dict['b_meta']['nsessions']))
            ses_markers = np.concatenate(np.array([[i+1]*ecounts for i in range(AQUA_dict['b_meta']['nsessions'])]))[:len(AQUA_dict['AQUA_vis'])]
            AQUA_dict['AQUA_vis']['session'] = ses_markers
            AQUA_dict['AQUA_rew']['session'] = ses_markers
        except KeyError: # if the animal is in the modelling data but not real, skip over. Won't be used anyway 
            continue 

        AQUA_dicts[data[i]['anim_id'][0]] = AQUA_dict

    AQUA_dicts = add_session_marker(AQUA_dicts,data_dict)

    return GLM_dfs, AQUA_dicts




prob_path = '/Users/grimal/Dropbox/100hours/hexaport/photometry_exp/reviewer_response/other_data/'

def import_GLMmat_prob(prob_path,data_dict,regions=['NAc']):
    '''
    Similar to the above function but specific to probability GLM data
    Doesn't produce AQUA_dict
    
    '''
    aGLM_dfs = []
    for region in regions:
        mat_data = loadmat(prob_path+'GLM_export_ProbabilityData_'+region)
        data     = mat_data['GLM_export'][0]
        
        GLM_dfs = []
        for i in range(len(data)):
            # variables for GLM analyses 
            GLM_dict               = {}
            GLM_dict['mouse']      = data[i]['anim_id'][0][:6]
            GLM_dict['alpha']      = np.concatenate(data[i]['alpha'])
            GLM_dict['AQUA_rpe']   = np.concatenate(data[i]['AQUA_rpe'])
            GLM_dict['port_id']    = np.concatenate(data[i]['port_id'])
            GLM_df                 = pd.DataFrame(GLM_dict)
            GLM_dfs.append(GLM_df)
            
        aGLM_dfs.append(GLM_dfs)
        
    return list(chain.from_iterable(aGLM_dfs))



### IMPORTING DATA SPECIFIC TO LAST FIGURE 

KLpath = '/Users/grimal/Dropbox (Personal)/100hours/hexaport/photometry_exp/other_data/last_fig'

def import_KL_AQUA_data(KLpath):
    '''
    File 1:
        where rows are animals and columns are da data 4 columns [end 1, start 2, end 2, start 3] and then 4 columns for KL means matched to these time periods.
    File 2:
        In file *2 there is a struct alpha_and_kl_6pg5_aqua where I have corresponding alpha and KL values from start 2->start3 
        but the dimensions are now rows = increasing scaling of alpha, columns are repetitions of the model simulation to see variance 
        in KL estimates based upon sampling (with a little added jitter in alpha because that seemed reasonable noise model)
    '''

    # first file: DA and KL 
    DA_KL = loadmat(path+'/ForLauraFinalFig')
    DA_KL = pd.DataFrame(DA_KL['da_and_kl_dat_nac'])
    
    # second file 
    alpha_KL = loadmat(path+'/ForLauraFinalFig2')['alpha_and_kl_6pg5_aqua'][0][0]
    alpha_KL = [pd.DataFrame(matrix) for matrix in alpha_KL] # alpha is first three (start 2, end 2, start 3) and then KL is last three 
    
    return DA_KL, alpha_KL



### IMPORTING DATA FOR MODEL COMPARISONS 


comp_file = '/Users/grimal/Dropbox (Personal)/100hours/hexaport/photometry_exp/other_data/SessionDataCompareModels.mat'

def import_modelcomp(comp_file):
    '''
    For specific plots in fig 6 to do with the AQUA model and comparison to other models. 
    There's more in mc, but not saved at the moment/don't know what it is
        
    '''
    mc   = loadmat(comp_file)['session']
    r2   = pd.DataFrame(mc[0][1][0][0][0][0]) # 16 rows, 9 columns
    inc  = pd.DataFrame(mc[0][1][0][0][0][1])# 16 rows, 9 columns
    tobsub    = pd.DataFrame(mc[0][1][0][0][0][5])
    delta_rew = pd.DataFrame(mc[0][1][0][0][0][3])-np.array(tobsub) # 16 rows, 9 columns
    
    cols = [mc[0][1][0][0][0][-3][0][i][0] for i in range(len(mc[0][1][0][0][0][-3][0]))]
    
    r2.columns        = cols
    inc.columns       = cols
    delta_rew.columns = cols
    
    # alpha values for optimal fit: [peak value, decay tau (in visits), dc offset]
    opt_fit   = pd.DataFrame(mc['all_coms'][0][1])
    opt_alpha = pd.concat([opt_fit[2] + (opt_fit[0]*np.exp(-ivis/opt_fit[1])) for ivis in range(0,1000)],axis=1)
    
    
    return r2, inc, delta_rew, opt_alpha 

### OTHER SUPPORT FUNCS

def add_session_marker(AQUA_dicts,data_dict):
    '''
    In original form, model data doesn't have session markers and is continuous.
    So here adding session number as a column to AQUA_dicts (can also add to aGLM_dfs in future using similar method)
    Also only adding to data that is in the range of visits, not time (time already has session markers)
    '''
    for mouse in data_dict.keys():
        mAQUA           = AQUA_dicts[mouse]
        mAQUA_rew_invis = mAQUA['AQUA_rew_invis']
        mAQUA_vis_invis = mAQUA['AQUA_vis_invis']
        
        # double check number of visits in behaviour and AQUA match 
        alick_df        = data_dict[mouse]['conc']['mlick_df']
        avisits         = alick_df[alick_df['unique_visit']==1].reset_index(drop=True)

        # add session column to AQUA
        mAQUA_rew_invis['session_number'] = avisits['session_number']
        mAQUA_vis_invis['session_number'] = avisits['session_number']
        mAQUA['AQUA_rew_invis']           = mAQUA_rew_invis
        mAQUA['AQUA_vis_invis']           = mAQUA_vis_invis
        
        mAQUA['AQUA_rew_invis']           = mAQUA['AQUA_rew_invis'].dropna()
        mAQUA['AQUA_vis_invis']           = mAQUA['AQUA_vis_invis'].dropna()
        if len(mAQUA['AQUA_rew_invis'])!=len(avisits):
            print('continued discrepancy')

    return AQUA_dicts 

        
def get_visits_AQUA(mAQUA_vis,lick_df,ses_n,bmeta):
    '''
    Get data from AQUA and reconfigure to standard visits format
    Ie one row per visit, with port and event_time and port_rank as columns
    In pd.dataframe
    AQUA data has to be for a specific mouse and should be just visits - not rewards or bmeta 
    Ie it should be a single dataframe that has columns for each port and rows for each visit, for one session
    '''
    beh_visits                = lick_df[lick_df['unique_visit']==1][['port','event_time']]
    AQUA_visits               = pd.DataFrame(mAQUA_vis[mAQUA_vis['session_number']==ses_n].iloc[:,0:6].idxmax(axis=1)+1).reset_index(drop=True)
    AQUA_visits['event_time'] = beh_visits['event_time'].reset_index(drop=True)
    AQUA_visits               = AQUA_visits.rename(columns={AQUA_visits.columns[0]: 'port'})
    
    port_dict                 = sf.port_2_rank(bmeta,ses_n)
    port_ranks               = np.array([int(port_dict[vis]) for vis in AQUA_visits['port']])
    AQUA_visits['port_rank'] = port_ranks
    
    return AQUA_visits        
        
            
def get_rews_AQUA(mAQUA_rews,lick_df,ses_n):
    '''
    Get data from AQUA and reconfigure to standard visits format
    Ie one row per visit, with port and event_time as columns
    In pd.dataframe
    AQUA data has to be for a specific mouse and should be just visits - not rewards or bmeta 
    Ie it should be a single dataframe that has columns for each port and rows for each visit, for one session
    '''
    beh_visits                = lick_df[lick_df['unique_visit']==1][['port','event_time']]
    AQUA_rews                 = pd.DataFrame(mAQUA_rews[mAQUA_rews['session_number']==ses_n].iloc[:,0:6].max(axis=1)).reset_index(drop=True)
    AQUA_rews['event_time']   = beh_visits['event_time'].reset_index(drop=True)
    AQUA_rews                 = AQUA_rews.rename(columns={AQUA_rews.columns[0]: 'rewarded'})

    return AQUA_rews

    
    
    
    
    
    
    
    
    
    
    