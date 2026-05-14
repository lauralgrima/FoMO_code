import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def opto_updaterule(opto_filepath,plot=True):     
    '''
    Imports .mat file for simulations of the update function only 
    '''
    
    full_data = scipy.io.loadmat(opto_filepath)['Summary_UpdateRule_OptoSim'][0][0]  
    cntl_rew_probs = pd.DataFrame(full_data[1])
    cntl_vis_probs  = pd.DataFrame(full_data[2])

    worst_rew_probs = pd.DataFrame(full_data[3])
    worst_vis_probs  = pd.DataFrame(full_data[4])
    
    sall_rew_probs = pd.DataFrame(full_data[5])
    sall_vis_probs  = pd.DataFrame(full_data[6])
    
    # take the ratio
    def new_matching(rew_probs,vis_probs): # just the same as previous matching function but configured for simulations
        all_match_vis, all_match_rew = [],[]
        for sim in range(len(rew_probs)):
            rew = rew_probs.iloc[sim]
            vis = vis_probs.iloc[sim]
            
            match_rew = np.array([rew[i]/(np.sum(rew)-rew[i]) for i in range(len(rew))])
            match_vis = np.array([vis[i]/(np.sum(vis)-vis[i]) for i in range(len(vis))])
            
            # remove zeros or nans
            match_vis = np.nan_to_num(match_vis,nan=0.001)
            match_vis[match_vis==float(0)]=0.001 
            match_rew = np.nan_to_num(match_rew,nan=0.001)
            match_rew[match_rew==float(0)]=0.001 
            
            all_match_vis.append(match_vis),all_match_rew.append(match_rew)
        return all_match_vis,all_match_rew
    
    cntl_match_vis,cntl_match_rew = new_matching(cntl_rew_probs,cntl_vis_probs)
    worst_match_vis,worst_match_rew = new_matching(worst_rew_probs,worst_vis_probs)
    sall_match_vis,sall_match_rew = new_matching(sall_rew_probs,sall_vis_probs)
    
    if plot:
    
        fig, ax = plt.subplots(1,1,figsize=(3.8,3.5))
        hfont = {'fontname':'Helvetica'}
        
        # for real data, -2.5 to 1.0
        ax.plot(np.arange(-1.5,0.6),np.arange(-1.5,0.6),'gray',alpha=0.3,linewidth=2,zorder=-20)
        for sim in range(len(cntl_rew_probs)):
    
            ax.scatter(np.log10(cntl_match_rew[sim]),np.log10(cntl_match_vis[sim]),color='lightgrey',s=15,zorder=-10)
            ax.scatter(np.log10(worst_match_rew[sim]),np.log10(worst_match_vis[sim]),color=['navy','navy','deepskyblue','deepskyblue','deepskyblue','deepskyblue'],s=15)
            
        # for real data, -2.5 to 0.5
        ax.set_xlim([-1.5,0])
        ax.set_ylim([-1.5,0])
        ax.set_xticks([-1.5,-1.0,-0.5,0])
        ax.set_yticks([-1.5,-1.0,-0.5,0])
        ax.set_yticklabels([-1.5,-1.0,-0.5,0])
        ax.set_xticklabels([-1.5,-1.0,-0.5,0])
    
        ax.set_xlabel('log$_{10}$ reward ratio',**hfont,fontsize='18')
        ax.set_ylabel('log$_{10}$ choice ratio',**hfont,fontsize='18')
        ax.spines[['right','top']].set_visible(False)
        ax.spines['left'].set_position(('outward', 10))
        ax.spines['bottom'].set_position(('outward', 10))
        ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
        ax.spines[['left','bottom']].set_linewidth(1.5)
            
        fig.tight_layout()
        plt.show()
