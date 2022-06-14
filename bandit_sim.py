import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

probabilities = [0.2, 0.3, 0.5, 0.5, 0.7, 0.8]
probabilities = [0.2, 0.8]
n_decisions = 100
policy = 'win-stay lose-shift'

def n_armed_bandit(probabilities,n_decisions,policy,plot=True,savefigs=False):
    
    options = list(range(0,len(probabilities)))
    
    rew_history        = []
    choice_history     = []
    prop_rew_0_all     = []
    prop_rew_1_all     = []
    for decision in range(0,n_decisions):
        
        if policy == 'random':
            choice   = random.choices(options,np.repeat(1/len(options),len(options)))[0]
            rewarded = random.choices([1,0],weights=(probabilities[choice],1-probabilities[choice]))[0]
            
        elif policy == 'win-stay lose-shift':
            if len(rew_history) == 0: # first decision, choose randomly 
                choice = random.choices(options,np.repeat(1/len(options),len(options)))[0]
            else:
                if rew_history[decision-1] == 1: # if previous choice was rewarded
                    choice = choice_history[decision-1] # repeat previous choice 
                else:
                    prev_choice = choice_history[decision-1]
                    remaining_options = [option for option in options if option != prev_choice]
                    choice = random.choices(remaining_options)[0]
            rewarded = random.choices([1,0],weights=(probabilities[choice],1-probabilities[choice]))[0]
            
        # elif polocy == 'matching':
        #     if len(rew_history) == 0: # first decision, choose randomly 
        #         choice = random.choices(options,np.repeat(1/len(options),len(options)))[0]
        #     else:
                
            
                    
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
        rew_history.append(rewarded)

        bandit_dict = {'choice':choice_history,'rewarded':rew_history}
        bandit_df = pd.DataFrame(bandit_dict) 
    
    cumulative_rew = np.cumsum(rew_history)
    
    if plot:
        
        fig,ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [4, 0.5]},figsize=(10,2))
        hfont = {'fontname':'Helvetica'}    
        fig.suptitle(policy,fontsize="18",**hfont,color='darkorange')   
        
        colors = list(bandit_df['rewarded'])
        scatter = ax[0].scatter(x=bandit_df.index,y=bandit_df['choice'],alpha=0.4,c=colors,cmap='bwr_r')
        ax[0].tick_params('both',length=0)
        ax[0].set_yticks(range(-1,bandit_df['choice'].max()+2))
        ax[0].set_yticklabels([])
        ax[0].set_xticks([0,100])
        ax[0].set_xticklabels([0,100],fontsize='12')
        ax[0].set_ylabel('chosen option',fontsize='14',**hfont)
        ax[0].set_xlabel('choice number',fontsize='14',**hfont)
        ax[0].legend(*scatter.legend_elements(),loc='lower right')
        for axis in ['top','bottom','left','right']:
            ax[0].spines[axis].set_linewidth(1) 
        
        ax[1].plot(bandit_df.index,cumulative_rew,color='blue')
        ax[1].tick_params('both',length=0)
        ax[1].set_ylim([0,100])
        ax[1].set_xticks([0,100])
        ax[1].set_yticks([0,100])
        ax[1].set_xticklabels([0,100],fontsize='12')
        ax[1].set_yticklabels([0,100],fontsize='12')
        ax[1].set_ylabel('cumulative rew.',fontsize='14',**hfont)
        ax[1].set_xlabel('choice number',fontsize='14',**hfont)
        for axis in ['top','bottom','left','right']:
            ax[1].spines[axis].set_linewidth(1) 

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])   
        
    if savefigs:
        plt.savefig(policy,dpi=300)


    return(bandit_df,cumulative_rew)
        

# plotting cumulative reward separately
    
# fig,ax = plt.subplots(1,1,figsize=(3,3))
# ax.tick_params('both',length=0)
# ax.set_ylim([0,100])
# ax.set_xticks([0,100])
# ax.set_yticks([0,100])
# ax.set_xticklabels([0,100],fontsize='12')
# ax.set_yticklabels([0,100],fontsize='12')
# ax.set_ylabel('cumulative rew.',fontsize='14',**hfont)
# ax.set_xlabel('choice number',fontsize='14',**hfont)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(1) 
# ax.plot(bandit_df.index,cum_random,color='darkorange',linewidth=2)
# fig.tight_layout(rect=[0, 0.03, 1, 0.95])   
# ax.plot(bandit_df.index,cum_wsls,color='green',linewidth=2)














