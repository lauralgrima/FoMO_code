import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

probabilities = [0.8, 0.2, 0.3]
n_decisions = 100
policy = 'random'

def n_armed_bandit(probabilities,n_decisions,policy,plot=True):
    
    options = list(range(0,len(probabilities)))
    
    rew_history        = []
    choice_history     = []
    rew_option_history = []
    prop_rew_0_all     = []
    prop_rew_1_all     = []
    for decision in range(0,n_decisions):
        
        # give which option is rewarded 
        rewarded_option = random.choices(options,weights=probabilities)[0]
        rew_option_history.append(rewarded_option)
        
        if policy == 'random':
            choice = random.choices(options,np.repeat(1/len(options),len(options)))[0]
            
        elif policy == 'wsls':
            if len(rew_option_history) == 1: # first decision, choose randomly 
                choice = random.choices(options,np.repeat(1/len(options),len(options)))[0]
            else:
                if rew_history[decision-1] == 1: # if previous choice was rewarded
                    choice = choice_history[decision-1] # repeat previous choice 
                else:
                    prev_choice = choice_history[decision-1]
                    remaining_options = [option for option in options if option != prev_choice]
                    choice = random.choices(remaining_options)
        # elif policy == 'matching':
        #     if len(rew_option_history) == 1: # first decision, choose randomly 
        #         choice = random.choices(options,weights=[0.5,0.5])[0]
        #     else: 
        #         # calculate proportion of choices for each option that have been rewarded 
        #         if len(bandit_df.loc[bandit_df['choice']==0]) == 0: # option hasn't been chosen yet 
        #             prop_rew_0 = 0.5
        #         else: 
        #             prop_rew_0 = (len((bandit_df.loc[(bandit_df['choice']==0) & (bandit_df['rewarded']==1)])))/(len(bandit_df.loc[bandit_df['choice']==0]))
        #         if len(bandit_df.loc[bandit_df['choice']==1]) == 0:
        #             prop_rew_1 = 0.5
        #         else: 
        #             prop_rew_1 = (len((bandit_df.loc[(bandit_df['choice']==1) & (bandit_df['rewarded']==1)])))/(len(bandit_df.loc[bandit_df['choice']==1]))
        #         prop_rew_0_all.append(prop_rew_0)
        #         prop_rew_1_all.append(prop_rew_1)
        #         if prop_rew_0 > prop_rew_1:
        #             choice = 0
        #         else:
        #             choice = 1
                    
        choice_history.append(choice)
        
        if choice == rewarded_option:
            reward = 1
        else:
            reward = 0
        rew_history.append(reward)
        
        bandit_dict = {'rewarded_option':rew_option_history,'choice':choice_history,'rewarded':rew_history}
        bandit_df = pd.DataFrame(bandit_dict) 
        print(bandit_df)
    
    cumulative_rew = np.cumsum(rew_history)

    
    if plot:
        plt.plot(cumulative_rew)
        plt.xlabel('decision number')
        plt.ylabel('cumulative reward')
        plt.title('win stay lose shift')
        
        
        
        

