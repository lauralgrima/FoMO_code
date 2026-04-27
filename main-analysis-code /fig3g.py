import matplotlib.pyplot as plt
import support_funcs as sf
import numpy as np
import pandas as pd
import seaborn as sns
import regression_variables_behavior as rv
import import_mat as im
from sklearn import linear_model
from scipy import stats
from sklearn.model_selection import cross_val_score, train_test_split


def MULTImultinomial_regr(data_dict,mat_GLM_path,ses_n=1,cv=5,pred_approach='train_test',plot=True):
    """
    Run multinomial logistic regression across mice for one session.
    
    Fits choice-prediction models for both experimental behavior and AQUA/model
    behavior, then collects regression coefficients, intercepts, prediction
    accuracy, and cross-validation scores across mice. Optionally plots summed
    absolute coefficients for selected predictors.
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral data.
    mat_GLM_path : str
        Path to GLM/AQUA model data.
    ses_n : int, default=1
        Session number to analyze.
    cv : int, default=5
        Number of cross-validation folds.
    pred_approach : str, default='train_test'
        Prediction approach passed to `multinomial_logr`.
    plot : bool, default=True
        Whether to plot coefficient summaries.
    """
    
    predictors      = ['n1_rew','n2_rew','n3_rew','n1_choice','n2_choice','n3_choice','n1_int','n2_int','n3_int','ratio_match','cumu_rew_port','time','dist']
    plot_predictors = ['n1_rew','n1_choice','n1_int','ratio_match','cumu_rew_port','time','dist']
    pred_labels     = ['prev. rew.','stay/leave','WSLS','matching','cumu. rew.','time','distance']

    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    subset_dict     = {**subset_dict_NAc,**subset_dict_DMS}
        
    _,AQUA_beh_NAc  = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS  = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_beh        = {**AQUA_beh_DMS,**AQUA_beh_NAc}
        
    acoef,aintercept,aperc_pred,acv_mean,acv_sem = [],[],[],[],[]
    aAQUAcoef,aAQUAintercept,aAQUAperc_pred,aAQUAcv_mean,aAQUAcv_sem = [],[],[],[],[]
    for mouse in subset_dict:
        print(mouse)
        
        # real data
        lick_df,_,vid_df,bmeta,_                       = sf.extract_data(data_dict,mouse,ses_n) 
        
        mouse_regr,perc_pred,cv_mean,cv_sem            = multinomial_logr(subset_dict,lick_df,vid_df,bmeta,cv,ses_n,predictors,data_type='real',pred_approach=pred_approach)
        AQUA_regr,AQUAperc_pred,AQUAcv_mean,AQUAcv_sem = multinomial_logr(subset_dict,[lick_df,AQUA_beh[mouse]],vid_df,bmeta,cv,ses_n,predictors,data_type='AQUA',pred_approach=pred_approach)
        if isinstance(mouse_regr,linear_model._logistic.LogisticRegression): 
            coef           = mouse_regr.coef_
            AQUA_coef      = AQUA_regr.coef_
            intercept      = mouse_regr.intercept_
            AQUA_intercept = AQUA_regr.intercept_
            pd_coef        = pd.DataFrame(coef,columns=predictors)
            pd_coef        = pd_coef[plot_predictors]
            
            acoef.append(pd_coef),aintercept.append(intercept),aperc_pred.append(perc_pred),acv_mean.append(cv_mean),acv_sem.append(cv_sem)
            aAQUAcoef.append(pd.DataFrame(AQUA_coef,columns=predictors)),aAQUAintercept.append(AQUA_intercept),aAQUAperc_pred.append(AQUAperc_pred),aAQUAcv_mean.append(AQUAcv_mean),aAQUAcv_sem.append(AQUAcv_sem)
            
    if plot: 
        plot_coef(acoef,aAQUAcoef,plot_predictors,pred_labels)
        plot_cv_accuracy(acv_mean,aAQUAcv_mean)


def multinomial_logr(data_dict,lick_df,vid_df,bmeta,cv,ses_n,predictors,data_type='real',pred_approach='train_test'):
    """
    Run multinomial logistic regression to predict port choices.
    
    Fits a logistic regression model using the specified predictors to predict
    observed port choices for real or AQUA data. Prediction accuracy is computed
    using either a train/test split or the full dataset, with optional
    cross-validation for real data.
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session data.
    lick_df : pandas.DataFrame or list
        Behavioral data. For AQUA, a list containing [lick_df, AQUA_beh_mouse].
    vid_df : pandas.DataFrame
        Video-derived predictor data.
    bmeta : object
        Behavioral metadata.
    cv : int
        Number of cross-validation folds; use 0 to skip.
    ses_n : int
        Session number.
    predictors : list
        Predictor names to include.
    data_type : {'real', 'AQUA'}, default='real'
        Type of data to analyze.
    pred_approach : str, default='train_test'
        If 'train_test', evaluates on held-out data; otherwise fits and predicts
        on the full dataset.
    
    Returns
    -------
    tuple
        Fitted model, prediction accuracy, cross-validation mean, and
        cross-validation SEM.
    """

    if data_type == 'AQUA':
       AQUA_beh_mouse = lick_df[1]
       lick_df        = lick_df[0]
   
    # observed choices 
    if data_type == 'real':
        obs_choices = rv.gen_observed(lick_df,bmeta)['port_ranks']
    elif data_type == 'AQUA':
        obs_choices = im.get_visits_AQUA(AQUA_beh_mouse['AQUA_vis_invis'],lick_df,ses_n,bmeta)['port_rank']

    if data_type == 'real':
        predictor_matrix = rv.gen_predictors(lick_df,bmeta,ses_n,predictors,data_type,scale=True)
    elif data_type == 'AQUA':
        predictor_matrix = rv.gen_predictors([lick_df,AQUA_beh_mouse],bmeta,ses_n,predictors,data_type,scale=True)
        if 'travel' in predictors:
            predictor_matrix = predictor_matrix.drop('travel',axis=1)

    # remove nan rows 
    inan             = np.unique(np.where(predictor_matrix.isna())[0])
    predictor_matrix = predictor_matrix[~predictor_matrix.index.isin(inan)].reset_index(drop=True)
    obs_choices      = obs_choices[~obs_choices.index.isin(inan)].reset_index(drop=True)

    # set up and run log regression model
    regr       = linear_model.LogisticRegression(max_iter=500)
    
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
    try:
        cv_mean,cv_sem = cross_val(regr,predictor_matrix,obs_choices,cv)
    except ValueError:
        cv_mean,cv_sem = [],[]
        mouse_regr     = []
            
    return mouse_regr,perc_pred,cv_mean,cv_sem


def cross_val(regr,predictor_matrix,observed_choices,cv):
    '''
    Cross validation of logistic regression - gives an accuracy score (i.e. how well the model predicts real choices).
    Variables:
        regr     :     scipy.linear_model.LogisticRegression model
        predictor_matrix    :    pd.DataFrame of predictors, one per column
        observed_choices    :    pd.Series of choices at a given port (1 or 0)
        cv                  :    number of cross validation runs. Set as integer. Usually 5 or 10
    Returns mean and SEM of all cv runs. 
    '''
    scores = cross_val_score(regr,predictor_matrix,observed_choices,cv=cv,scoring='accuracy')
    return np.mean(scores), stats.sem(scores)


### PLOTTING

def plot_coef(coef_dfs,AQUAcoef_dfs,plot_predictors,pred_labels):
    """
    Plot summed absolute regression coefficients across predictors.
    
    For each predictor, computes the sum of absolute coefficients across all
    ports (1–6), averaged across datasets, and plots these values with SEM
    error bars for both experimental and AQUA data.
    
    Parameters
    ----------
    coef_dfs : list of pandas.DataFrame
        Coefficient DataFrames for experimental data.
    AQUAcoef_dfs : list of pandas.DataFrame
        Coefficient DataFrames for AQUA/model data.
    plot_predictors : list
        Subset of predictors to include in the plot.
    pred_labels : list
        Labels corresponding to `plot_predictors` for x-axis display.
    """
    
    # plotting summed absolute values across ports
    fig,ax    = plt.subplots(figsize=(3.8,4))
    hfont     = {'fontname':'Arial'}
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    colors = ['black','lightseagreen']
    
    for i,dfs in enumerate([coef_dfs,AQUAcoef_dfs]):

        all_coefs = pd.concat(dfs)
        all_coefs = all_coefs[plot_predictors]

        amean_coefs = pd.concat([all_coefs[all_coefs.index==j].mean() for j in range(0,6)], axis=1)
        asem_coefs  = pd.concat([all_coefs[all_coefs.index==j].sem()  for j in range(0,6)], axis=1)
    
        abs_means = amean_coefs.abs().sum(axis=1)
        abs_sems  = asem_coefs.abs().sum(axis=1)
        ax.spines[['left','bottom']].set_position(('outward', 10))
        ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
        ax.scatter(x=range(0,len(abs_means)),y=abs_means,color=colors[i])
        ax.errorbar(x=range(0,len(abs_means)),y=abs_means,yerr=abs_sems,fmt='none',color=colors[i])

    ax.set_ylabel(r'abs. sum $\beta$' + '\ncoefficient',**hfont,fontsize=18)
    ax.set_xticks(range(0,len(abs_means)))
    ax.set_xticklabels(pred_labels,rotation=45,ha='right',**hfont,fontsize=18)
    ax.set_yticks([0,10,20])
    ax.set_yticklabels([0,10,20])
    ax.spines[['right','top']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.xaxis.set_tick_params(width=1.5)

    fig.tight_layout()
    
    
def plot_cv_accuracy(cv_means1,cv_means2):
    """
    Plot cross-validated prediction accuracy for two sessions.
    
    Shows paired distributions of cross-validation mean accuracies, with individual
    points and lines connecting matched mice across sessions.
    
    Parameters
    ----------
    cv_means1 : array-like
        Cross-validation mean accuracies for the first session.
    cv_means2 : array-like
        Cross-validation mean accuracies for the second session.
    """
    
    hfont = {'fontname':'Arial'} 
    fig,ax = plt.subplots(figsize=(3,3.5))

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    sns.boxplot(y=cv_means1,x=0,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False,zorder=10)
    sns.stripplot(y=cv_means1,x=0,edgecolor='gainsboro',facecolor='white',linewidth=2)
    
    sns.boxplot(y=cv_means2,x=0.2,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False,zorder=10)
    sns.stripplot(y=cv_means2,x=0.2,edgecolor='gainsboro',facecolor='white',linewidth=2)
    
    [sns.lineplot(x=[0,1],y=[cv_means1[i],cv_means2[i]],ax=ax,color='gainsboro',linewidth=2,zorder=-10) for i in range(0,len(cv_means1))]
    
    ax.set_ylabel('C.V. predicted choices (%)',**hfont,fontsize=18)
    
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_yticklabels([0,20,40,60,80,100],fontsize=18)
    ax.set_xlim([-1,5])
    ax.set_ylim([0,1.0])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)

    fig.tight_layout()

