import numpy as np
import matplotlib.pyplot as plt
import preprocessing.support_funcs as sf
import pandas as pd
from scipy.stats import spearmanr


def MULTI_RW_DA(data_dict,predictors,ses_n,tw_start=-0.5,tw_length=1.5,epsilon=0.1,plot=True):
    """
    Compute RW-DA correlations across mice.
    
    Runs `RW_DA` for each mouse to correlate dopamine response magnitude with
    Rescorla-Wagner RPE and value. Can analyze one session or average correlations
    across all sessions per mouse.
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral and photometry data.
    predictors : list
        Currently unused.
    ses_n : int or str
        Session number to analyze, or 'all' to average across all sessions.
    tw_start : float, default=-0.5
        Window start time relative to visit.
    tw_length : float, default=1.5
        Window length in seconds.
    epsilon : float, default=0.1
        Learning rate for the Rescorla-Wagner model.
    plot : bool, default=True
        Currently unused.
    
    Returns
    -------
    None
    """
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
            
def RW_DA(lick_df,photo_df,pmeta,ses_n,by_port=False,tw_start=-0.5,tw_length=1.5,epsilon=0.1,plot=False):
    """
    Correlate Rescorla-Wagner variables with dopamine response magnitude.
    
    Computes Rescorla-Wagner value and prediction error for each port, extracts
    visit-aligned photometry responses, summarizes each response by AUC, and
    correlates normalized DA responses with RW value or RPE.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing unique visits, ports, rewards, and photo indices.
    photo_df : pandas.DataFrame
        Photometry data containing a 'signal' column.
    pmeta : dict-like
        Photometry metadata containing sampling rates.
    ses_n : int
        Session number.
    by_port : bool, default=False
        Whether to compute correlations separately by port.
    tw_start : float, default=-0.5
        Window start time relative to visit.
    tw_length : float, default=1.5
        Window length in seconds.
    epsilon : float, default=0.1
        Learning rate for the Rescorla-Wagner model.
    plot : bool, default=False
        Whether to plot RW variables and DA responses.
    
    Returns
    -------
    tuple
        Spearman correlations for RPE and value versus DA response.
    """
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
    """
    Compute Rescorla-Wagner value and prediction error by port.
    
    Updates reward expectation independently for each port across visits, using a
    fixed learning rate. Values and prediction errors are then recombined in the
    original visit order.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing 'unique_visit', 'port', and 'rewarded'.
    epsilon : float
        Learning rate for the Rescorla-Wagner update.
    plot : bool, default=False
        Whether to plot value and prediction error traces for each port.
    
    Returns
    -------
    pandas.DataFrame
        Visit-level Rescorla-Wagner outputs with columns:
        'vals', 'rew', 'pes', and 'vcount'.
    """

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
