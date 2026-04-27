import os
import matplotlib.pyplot as plt
import support_funcs as sf
import numpy as np
import pandas as pd
import seaborn as sns
import fig1g_ss_diverge as div
import regression_variables as rv
import import_mat as im
from sklearn import linear_model
from scipy import stats
from sklearn.model_selection import cross_val_score, train_test_split
from matplotlib.pyplot import cm






### MODEL COMPARISONS


# full model values TRYING RATIO MATCH
predictors      = ['n1_rew','n2_rew','n3_rew','n1_choice','n2_choice','n3_choice','n1_int','n2_int','n3_int','ratio_match','cumu_rew_port','time','dist']#,'travel']#,'rpe','value']
plot_predictors = ['n1_rew','n1_choice','n1_int','ratio_match','cumu_rew_port','time','dist','travel']
pred_labels     = ['n-1 rew.','n-1 choice','n-1 int.','matching','cumu. rew','time','distance','travel']

# FOR AQUA

predictors      = ['n1_rew','n2_rew','n3_rew','n1_choice','n2_choice','n3_choice','n1_int','n2_int','n3_int','ratio_match','cumu_rew_port','time','dist']
plot_predictors = ['n1_rew','n1_choice','n1_int','ratio_match','cumu_rew_port','time','dist']
pred_labels     = ['n-1 rew.','n-1 choice','n-1 int.','matching','cumu. rew','time','distance']




gm_predictors = ['n1_rew','n2_rew','n3_rew','n1_choice','n2_choice','n3_choice','n1_int','n2_int','n3_int','ratio_match','cumu_rew','time','dist']#,'rpe','value']


# for comp_models and plotting 
model_list  = [predictors,drop_matching,drop_cumu,drop_wsls,drop_shorthist,drop_time,drop_dist,drop_travel]
model_names = ['full','matching','cumu. rew','WSLS','short-term hist.','time','distance','travel']








def multinomial_logr(data_dict,lick_df,vid_df,bmeta,rank,cv,ses,data_type='real',pred_approach='train_test'):
    '''
    Runs a multinomial regression with predictor as all choices. 
    FOR REAL VS MODEL SIM DATA:
        Set data_type to 'real' to run on real mouse data.
        Set to 'sim' to run on simulated data. If doing this, lick_df should be the already reduced three column df (ports, rewarded, event_time)
        rather than a full dataframe. I.e. simdata_run.df 
        Make sure not to include travel as a predictor, otherwise all data will be removed as travel column is just NaN
    
    FOR AQUA:
        lick_df should be a list of two items: lick_df, and AQUA data for a mouse 
        Function splits these up 
    
    Set CV to 0 to not do any cross validation. 
    Set rank to True to do this for port ranks, not physical ports per se. 
    Set pred_approach to train_test to split data into train_test datasets. Otherwise will train on the full dataset 

    Returns:
        coefficients       :   a numpy array of length 6, one for each port, each of length number of predictros 
        intercepts         :   a list of six values, one for each port
        perc_pred          :   percentage of correctly predicted choices
        cv_mean            :   mean cross-validated predicted choice
        cv_sem             :   SEM of cross-validated predicted choices 
    '''
    
    predictors      = ['n1_rew','n2_rew','n3_rew','n1_choice','n2_choice','n3_choice','n1_int','n2_int','n3_int','ratio_match','cumu_rew_port','time','dist']
    plot_predictors = ['n1_rew','n1_choice','n1_int','ratio_match','cumu_rew_port','time','dist']
    pred_labels     = ['n-1 rew.','n-1 choice','n-1 int.','matching','cumu. rew','time','distance']


    
    
    
    if data_type == 'AQUA':
       AQUA_beh_mouse = lick_df[1]
       lick_df        = lick_df[0]
   
    # observed choices 
    if data_type == 'real':
        obs_choices = rv.gen_observed(lick_df,bmeta)['port_ranks']
    elif data_type == 'sim':
        obs_choices = lick_df['port']
    elif data_type == 'AQUA':
        obs_choices = im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses,bmeta)['port_rank']

    if data_type == 'real':
        predictor_matrix = rv.gen_predictors(data_dict,lick_df,vid_df,bmeta,predictors,[ses],data_type,import_travels=True)
    elif data_type == 'AQUA':
        predictor_matrix = rv.gen_predictors(data_dict,[lick_df,AQUA_beh_mouse],vid_df,bmeta,predictors,[ses],data_type='AQUA',import_travels=True)
        if 'travel' in predictors:
            predictor_matrix = predictor_matrix.drop('travel',axis=1)
        

    # remove nan rows 
    inan             = np.unique(np.where(predictor_matrix.isna())[0])
    predictor_matrix = predictor_matrix[~predictor_matrix.index.isin(inan)].reset_index(drop=True)
    obs_choices      = obs_choices[~obs_choices.index.isin(inan)].reset_index(drop=True)

    # set up and run log regression model
    regr       = linear_model.LogisticRegression(multi_class='auto',max_iter=500) # auto means it will select ovr for binary data and multinomial for non-binary data 
    
    if pred_approach == 'train_test':
        # Split the dataset into training and testing sets
        try:
            X_train, X_test, y_train, y_test = train_test_split(predictor_matrix, obs_choices, test_size=0.2)
            mouse_regr                       = regr.fit(X_train,y_train)
            predicted_choices                = regr.predict(X_test)
        except ValueError:
            predicted_choices = []
    else:
        # train on whole dataset
        mouse_regr = regr.fit(predictor_matrix,obs_choices)
        predicted_choices = regr.predict(predictor_matrix)

    if len(predicted_choices)>0:
        if pred_approach == 'train_test':
            perc_pred = (predicted_choices==np.array(y_test)).astype(int).sum()/len(predicted_choices)
        else:
            perc_pred = (predicted_choices==np.array(obs_choices)).astype(int).sum()/len(predicted_choices)
    else: 
        perc_pred,mouse_regr = [],[]
    
    #cross validation
    if data_type == 'real':
        try:
            cv_mean,cv_sem = cross_val(regr,predictor_matrix,obs_choices,cv)
        except ValueError:
            cv_mean,cv_sem = [],[]
            mouse_regr     = []
    else:
        cv_mean,cv_sem = [],[]
            
    return mouse_regr,perc_pred,cv_mean,cv_sem



def plot_coef(coef_dfs,AQUAcoef_dfs,plot_predictors,pred_labels,comp,savefigs):
    '''
    For each predictor, plot:
        1. Absolute sum of normed coefficients (across all 6 ports)
    '''
    
    # plotting summed absolute values across ports
    fig,ax    = plt.subplots(figsize=(3.8,4))
    hfont     = {'fontname':'Arial'}
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    colors = ['lightseagreen','black']
    
    for i,dfs in enumerate([coef_dfs,AQUAcoef_dfs]):
    
        if not comp:
            all_coefs = pd.concat(dfs)
            all_coefs = all_coefs[plot_predictors]
        else:
            all_coefs = dfs
    
        amean_coefs = pd.concat([all_coefs[all_coefs.index==i].mean() for i in range(0,6)],axis=1) # average for each port
        asem_coefs  = pd.concat([all_coefs[all_coefs.index==i].sem() for i in range(0,6)],axis=1)
    
        abs_means = amean_coefs.abs().sum(axis=1)
        abs_sems  = asem_coefs.abs().sum(axis=1)
        ax.spines[['left','bottom']].set_position(('outward', 10))
        ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
        ax.scatter(x=range(0,len(abs_means)),y=abs_means,color=colors[i])
        ax.errorbar(x=range(0,len(abs_means)),y=abs_means,yerr=abs_sems,fmt='none',color=colors[i])
    if comp:
        ax.set_ylabel(r'$\Delta$ $\beta$ coefficient',**hfont,fontsize=18)
    else:
        ax.set_ylabel(r'abs. sum $\beta$' + '\ncoefficient',**hfont,fontsize=18)
    ax.set_xticks(range(0,len(abs_means)))
    ax.set_xticklabels(pred_labels,rotation=45,ha='right',**hfont,fontsize=18)
    if not comp:
        ax.set_yticks([0,10,20])
        ax.set_yticklabels([0,10,20])
    else:
        ax.set_yticks([-2,0,10])
        ax.set_yticklabels([-2,0,10])
    ax.spines[['right','top']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.xaxis.set_tick_params(width=1.5)

    fig.tight_layout()
    
    if savefigs:
        fig.savefig(os.path.join(savefigs,'fig6_AQUAcoef_fullval'),dpi=300,transparent=True)
        
    # plotting a single predictor split by port across x axis (mostly for showing the matching)
    fig2,ax2 = plt.subplots(figsize=(4,3.5))
    color    = cm.plasma(np.linspace(0.1,0.9,num=6)) # colour each port 
    
    predictor = 'ratio_match'
    
    ax2.axhline(y=0,color='lightgrey',zorder=-1)
    for port in range(0,6):
        sns.violinplot(y=all_coefs[all_coefs.index==port][predictor],x=port-0.4,ax=ax2,color=color[port],fill=False,width=0.4,split=True,inner=None,native_scale=True)
        sns.boxplot(y=all_coefs[all_coefs.index==port][predictor],x=port+0.1,ax=ax2,color=color[port],fill=False,width=0.2,showfliers=False,showcaps=False,native_scale=True)

    ax2.set_ylabel(r'mean $\beta$ coefficient',**hfont,fontsize=18)
    ax2.set_xticks([])
    ax2.set_xticklabels([])
    ax2.set_yticks([18,0,-10])
    ax2.set_ylim([-10,18])
    
    ax2.spines[['top','right','bottom']].set_visible(False)
    ax2.spines[['left']].set_linewidth(1.5)
    ax2.yaxis.set_tick_params(width=1.5)
    ax2.tick_params(axis='y',which='both', labelsize=18,direction='in')

    ax2.spines[['left']].set_position(('outward', 10))
        
    fig2.tight_layout()
    
    if savefigs:
        fig2.savefig(os.path.join(savefigs,'mean_beta_coef_match'),dpi=300,transparent=True)

