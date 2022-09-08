from scipy.special import softmax
import numpy as np
import pandas as pd
import model_data_df as mdf
import random

# import multi_lick_df file, 6PG5_NAc_conc.csv
multi_lick_df = pd.read_csv('/Users/grimal/Documents/GitHub/hexaport_model/eg_data/6PG5_NAc_conc.csv').iloc[:,1:]
video_df      = pd.read_csv('/Users/grimal/Documents/GitHub/hexaport_model/eg_data/video_df_6PG5_NAc_conc_2021-09-09-145129.csv').iloc[:,1:]

video_df1 = mdf.model_df(multi_lick_df,video_df,1) # huge df of all the relevant things 
policy_type = 'softmax' # softmax, greedy, e-greedy
belief_type = 'win_stay'
    

def model_choices(video_df1,policy_type,belief_type):
    
    max_rew_avail = sum(video_df1.iloc[:,8:13].sum())
    
    # TO ADD: start new kernel for randomness

    max_tsteps = len(video_df1)
    
    # creating dataframes for model output 
    hexa_model_visits  = pd.DataFrame(0,index=np.arange(len(video_df1)),columns=[1,2,3,4,5,6])
    hexa_model_rewards = pd.DataFrame(0,index=np.arange(len(video_df1)),columns=[1,2,3,4,5,6])
    p_reward           = pd.DataFrame(0,index=np.arange(len(video_df1)),columns=[1,2,3,4,5,6])
    reward_available   = pd.DataFrame(0,index=np.arange(len(video_df1)),columns=[1,2,3,4,5,6])
    
    checks             = video_df1.iloc[:,14:20]
    sample_logic       = checks.sum(axis=1)
    
    # maybe shorten to only 
    
    # run by videoframe timesteps   
    for t in range(1,max_tsteps):
        
        t = int(t)
        if sample_logic[t] == 1:
    
            
            
            
            
  #  for t in range(1,max_tsteps):
     #   p_reward.iloc[t,:] = p_reward.iloc[t-1,:] # keep track of reward probabilities across ports 
     #   reward_available.iloc[t,:] = video_df1.iloc[t,8:13]
     #   reward_available.iloc[t+1,:] = reward_available.iloc[t,:] # carry availability over to next check 
        
      #  for index, row in video_df1.iterrows():
        
        # should we check any port at this timepoint 
   #     if (checks.iloc[t] == 1).any():
            
            if policy_type == 'softmax':
                p_choice     = softmax((p_reward.iloc[t,:]/sum(p_reward.iloc[t,:])))
                checked_port = int([np.random.rand(1) > choice for choice in np.cumsum(p_choice)][-1]+1)
                
            elif policy_type == 'greedy':
                checked_port = max(p_reward.iloc[t,:])
                
            elif policy_type == 'e-greedy':
                if np.random.rand(1)>0.2:
                    checked_port = max(p_reward.iloc[t,:])
                else: 
                    checked_port = random.choice([1,2,3,4,5,6])
                
            elif policy_type == 'random':
                checked_port = random.choice([1,2,3,4,5,6])
                
            hexa_model_visits.iloc[t,checked_port]=1
            
            # was the check rewarded?
            if reward_available.iloc[t,checked_port]==1:
                hexa_model_rewards[t,checked_port]=1
                reward_available.iloc[t+1,:] = reward_available.iloc[t,:]
                reward_available.iloc[t+1,checked_port]=0
                yes_reward = 1
            else:
                yes_reward = 0
                
            # update belief (Pr(R|port,t)) according to different models 
                
            if belief_type == 'win_stay': # biased towards staying at current port after reward; visit with no reward explores
                p_reward.iloc[t,:] = 0.02
                if yes_reward:
                    p_reward[t,checked_port] = 0.9
                else: 
                    p_reward[t,checked_port] = 0.02
    
    print('Model rewards collected %i' %sum(hexa_model_rewards))
    print('Mouse rewards collected %i' %len(video_df1.loc[video_df1['rewarded']==1]))
    print('Max rewards available %i' %max_rew_avail)
    

    
    # check for nonzero rows
    p_reward.loc[~(p_reward==0).all(axis=1)]

    
    
    # matrix of port distances
    distances = [{'1':[0,14,18,70,72.2,65.5], '2':[14,0,22.8,56,65.5,42],'3':[18,22.8,0,72.2,70,56],'4':[70,56,72.2,0,18,22.8],'5':[72.2,65.5,70,18,0,14],'6':[65.5,42,56,22.8,14,0]}]
    distances_df = pd.DataFrame(distances)
    distances_df = distances_df.apply(pd.Series.explode)
    distances_df.index = [1,2,3,4,5,6]
    
