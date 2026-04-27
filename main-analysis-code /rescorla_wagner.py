import numpy as np
import matplotlib.pyplot as plt
import preprocessing.support_funcs as sf
import pandas as pd
from scipy.stats import spearmanr

# rescorla wagner calculations

### ACROSS MULTIPLE MICE 
def multi_RW_DA(data_dict,predictors,ses_n,tw_start=-0.5,tw_length=1.5,epsilon=0.1,plot=True,savefigs='/Users/grimal/Dropbox (Personal)/100hours/hexaport/photometry_exp/plots/initial_manu_plots'):
    '''
    Running RW_DA across mice to calculate correlation between DA peaks and either RPE or value from RW model.
    Set ses_n to string 'all' to calculate across all sessions. 
    '''
    all_rpe_corr,all_val_corr = [],[]
    for i,mouse in enumerate(list(data_dict.keys())):
        print(mouse)
        # if just for one session vs. all sessions 
        if isinstance(ses_n,int): 
            tot_ses = data_dict[mouse]['b_meta']['nsessions'] # total number of sessions for given mouse 
            if tot_ses>=ses_n: # if session is in range for given mouse 
                lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,ses_n) 
                rpe_corr,val_corr = RW_DA(lick_df,photo_df,pmeta,ses_n,by_port=False,tw_start=tw_start,tw_length=tw_length,epsilon=epsilon,plot=False)
                rpe_corr = rpe_corr[0]
                val_corr = val_corr[0]
        elif isinstance(ses_n,str): # run through all sessions for a mouse
            tot_ses = data_dict[mouse]['b_meta']['nsessions']
            arpe_corr,aval_corr = [],[]
            for i in range(1,tot_ses+1):
                lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,i) 
                rpe_corr,val_corr = RW_DA(lick_df,photo_df,pmeta,i,by_port=False,tw_start=tw_start,tw_length=tw_length,epsilon=epsilon,plot=False)
                arpe_corr.append(rpe_corr), aval_corr.append(val_corr) # animal
            rpe_corr = np.mean([arpe_corr[i][0] for i in range(0,len(arpe_corr))]) # average across sessions
            val_corr = np.mean([aval_corr[i][0] for i in range(0,len(aval_corr))])
        all_rpe_corr.append(rpe_corr)
        all_val_corr.append(val_corr)
            

### CALCULATIONS

# below function might be useful for plotting quartiles/response comparisons based on e.g. size of RPE 

def RW_DA(lick_df,photo_df,pmeta,ses_n,by_port=False,tw_start=-0.5,tw_length=1.5,epsilon=0.1,plot=False):
    '''
    Correlate RW RPE or value with DA AUC (normalised).
    Set by_port to True to do one correlation per port. Otherwise correlates across all ports. 
    '''
    values,rpes = rescorla_wagner(lick_df,epsilon,plot=False)
    photo       = photo_df['signal'].to_numpy()

    rpe_corr,val_corr = [],[]
    all_norm_neg_auc,all_port_rpes,all_norm_auc,all_port_vals = [],[],[],[]
    for iport in range(0,6):
        port_values = values[iport]
        port_rpes   = rpes[iport]
        
        # get DA responses - peaks or troughs - to rewarded/unrewarded visits 
        iphoto                = lick_df[(lick_df['unique_visit']==1)&(lick_df['port']==iport+1)]['photo_i']  
        photo_win,idx_ignored = sf.extract_window(iphoto,photo,pmeta['sampling_rate'][ses_n-1],tw_start=tw_start,tw_length=tw_length)
        sig_win               = np.vstack([photo[photo_win[x][0]:photo_win[x][1]] for x in range(len(photo_win))])
        sig_win0              = np.array([sig-sig[0] for sig in sig_win]) # zero 
        auc_sig               = np.array([np.trapz(sig,dx=5) for sig in sig_win0])
        norm_neg_auc          = sf.NormalizeNeg(auc_sig) # between -1 and 1 for rpe
        norm_auc              = sf.NormalizeData(auc_sig) # between 0 and 1 for value
        
        # use idx_ignored to drop values or RPEs of photometry trials that aren't being included
        if idx_ignored:
            port_values = np.delete(port_values,idx_ignored)
            port_rpes   = np.delete(port_rpes,idx_ignored)

        if plot:
            # plot RPE and values in order with DA responses 
            plt.figure(),plt.plot(port_rpes['pes']),plt.scatter(x=range(0,len(port_rpes)),y=norm_neg_auc),plt.xlabel('visits'),plt.ylabel('norm response'),plt.title('RPE')
            plt.figure(),plt.plot(port_values['vals']),plt.scatter(x=range(0,len(port_rpes)),y=norm_auc),plt.xlabel('visits'),plt.ylabel('norm response'),plt.title('value')

            # RPE: re-order for visualisation
            by_rpe      = sorted(zip(port_rpes['pes'],norm_neg_auc))
            sorted_rpes = np.array([by_rpe[i][0] for i in range(0,len(by_rpe))])
            sorted_sig  = [by_rpe[i][1] for i in range(0,len(by_rpe))]
            plt.figure(), plt.title('RPE x DA AUC response'), plt.xlabel('visit'),plt.ylabel('norm response')
            plt.scatter(range(0,len(sorted_rpes)),sorted_rpes,color='green',s=5)
            plt.scatter(range(0,len(sorted_sig)),sorted_sig,color='red',s=5)
            
            # value: re-order for visualisation
            by_val      = sorted(zip(port_values,norm_auc))
            sorted_vals = np.array([by_val[i][0] for i in range(0,len(by_val))])
            sorted_sig  = np.array([by_val[i][1] for i in range(0,len(by_val))])
            plt.figure(), plt.title('val x DA AUC response'), plt.xlabel('visit'),plt.ylabel('norm response')
            plt.scatter(range(0,len(sorted_vals)),sorted_vals,color='green',s=5)
            plt.scatter(range(0,len(sorted_sig)),sorted_sig,color='red',s=5)
            
            # matching correlation - add global?
            # port_vis = lick_df[(lick_df['unique_visit']==1)&(lick_df['port']==iport+1)]
            # prop_rew_lo = port_vis.cumsum()['rewarded']/list(range(0,len(port_vis)))
            # prop_rew_lo  = prop_rew_lo.replace(np.inf,np.nan)
            # by_match = sorted(zip(prop_rew_lo,norm_neg_auc))
            # sorted_match = np.array([by_match[i][0] for i in range(0,len(by_match))])
            # sorted_sig   = np.array([by_match[i][1] for i in range(0,len(by_match))])
            # plt.figure(), plt.title('match x DA AUC response'), plt.xlabel('visit'),plt.ylabel('norm response')
            # plt.scatter(range(0,len(sorted_match)),sorted_match,color='green',s=5)
            # plt.scatter(range(0,len(sorted_sig)),sorted_sig,color='red',s=5)


            
   #         prop_rew_lo  = obs.groupby('port').cumsum()['rewarded']/obs.groupby('port').cumcount() # proportion of checks at that port that have been rewarded
  #          prop_rew_glo = obs.groupby('port').cumsum()['rewarded']/obs['rewarded'].cumsum()
            
        if by_port:
            rpe_corr.append(spearmanr(norm_neg_auc,port_rpes))
            val_corr.append(spearmanr(norm_auc,port_values))
        else:
            all_norm_neg_auc.append(norm_neg_auc)
            all_port_rpes.append(port_rpes)
            all_norm_auc.append(norm_auc)
            all_port_vals.append(port_values)
    
    rpe_corr = spearmanr(np.concatenate(all_norm_neg_auc),np.concatenate(all_port_rpes))
    val_corr = spearmanr(np.concatenate(all_norm_auc),np.concatenate(all_port_vals))
    
    return rpe_corr,val_corr 
            

def rescorla_wagner(lick_df,epsilon,plot=False):
    '''
    Calculating RPE based on RW model, independent for each port. 
    Time is in the space of visits, not time per se. 
    Epsilon is learning rate, not fitted here. 
    Set plot to True to plot value and RPE for each port across visits. 

    Returns a dataframe with the following columns:
        vals  : value 
        rew   : whether visit was rewarded
        pes   : prediction errors
    This dataframe is integrated over all ports, i.e. is length of total number of visits and gives 

    '''

    # adding a visit count to lick_df, for rearranging values when combining again later
    visits           = lick_df[lick_df['unique_visit']==1].reset_index(drop=True)
    visits['vcount'] = np.arange(len(visits))
    
    RWs = []
    for port in range(1,7):
        # keeping track
        port_vcount = visits[visits['port']==port]['vcount'].reset_index(drop=True)
        
        # for a port, calculate reward expectation at a given visit 
        rews  = lick_df[(lick_df['unique_visit']==1)&(lick_df['port']==port)]['rewarded']
        rews  = rews.reset_index(drop=True)

        v = 1 # prediction - start at 1 as ports are baited 
        w = 1 # weight
        epsi = epsilon # learning rate
        delta = 0 # prediction error
        time = np.arange(0,len(rews),1) # time in visit space
        
        preds,pes = [],[]
        for i in range(0,time[-1]+1):
            r = rews.iloc[i]
            v = w
            delta = r-v # pe
            w = w+epsi*delta # update weight using RW
            preds.append(v) # value
            pes.append(delta) # pe
            
        if plot: 
            # plot expectation for a port
            plt.figure(),plt.plot(preds),plt.ylabel('estimated value'),plt.xlabel('visit number'),plt.title('est. value by visit - port %i' %port)
            # plot RPE
            plt.figure(),plt.plot(pes),plt.ylabel('RPE'),plt.xlabel('visit number'),plt.title('RW-calculated RPE, port %i' %port)
    
        RW            = pd.DataFrame(preds,columns=['vals'])
        RW['rew']     = rews
        RW['pes']     = pes
        RW['vcount']  = port_vcount

        RWs.append(RW)        

    return pd.concat(RWs).sort_values(by='vcount')
