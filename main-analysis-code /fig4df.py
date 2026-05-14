import support_funcs as sf
import regression_variables_photometry as rvp
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.preprocessing import StandardScaler
from pypalettes import load_cmap
from statsmodels.stats.multitest import multipletests


def MULTIplinr_whole(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,data_dict,ses_n=[1,2],tw_start=-2,tw_length=5,plot=True,rew_only=False,tobepred='full'):
    """
    Run photometry regression across NAc mice for selected sessions.
    
    For each mouse/session, fits a ridge regression predicting photometry responses
    from behavioral predictors, collects coefficients, R² values, and predictor
    performance drops, and optionally plots time-resolved coefficients or
    performance-drop summaries.
    
    Parameters
    ----------
    travels_filepath : str
        Path to saved travel-related data.
    div_from_rand_filepath : str
        Path to divergence-from-random data.
    mat_GLM_path : str
        Path to GLM/AQUA model data.
    GLM_prob_filepath : str
        Path to GLM probability data.
    data_dict : dict
        Mouse/session behavioral and photometry data.
    ses_n : list, default=[1, 2]
        Sessions to include.
    tw_start : float, default=-2
        Start of photometry window relative to visit.
    tw_length : float, default=5
        Length of photometry window.
    plot : bool, default=True
        Whether to plot results.
    rew_only : bool, default=False
        Whether to include only rewarded visits.
    tobepred : {'full', 'peak'}, default='full'
        Photometry target to predict.
    
    Returns
    -------
    None
    """
    subset_dict = sf.subset_mice(data_dict,config=1,region='NAc')
    
    if tobepred == 'full' and rew_only == False: # fig4d plots 
        photo_predictors = ['rewarded','port_quality','prev_port_qual','n1_rew','n1_choice','n1_int','after_rand','time','dist']
        beh_predictors   = photo_predictors 
        plot_predictors  = ['rewarded','port_quality','prev_port_qual','n1_choice','after_rand']
        predic_labels    = plot_predictors

    elif tobepred == 'full' and rew_only == True: # fig4f, top row plots 
        photo_predictors = ['port_quality','prev_port_qual','n1_choice','time','dist','alpha','AQUA_rpe']
        beh_predictors   = photo_predictors 
        plot_predictors  = ['alpha','AQUA_rpe']
        predic_labels    = plot_predictors 
        
    elif tobepred == 'peak' and rew_only == True: # fig4f, bottom row plots
        photo_predictors = ['port_quality','prev_port_qual','n1_choice','after_rand','time','dist','alpha','AQUA_rpe']
        beh_predictors   = photo_predictors 
        plot_predictors  = photo_predictors
        predic_labels    = photo_predictors
        
    nac_pcoefs1, nac_pcoefs2 = [],[]
    r2s1, r2s2 = [],[]
    perf_drops1,perf_drops2 = [],[]
    for mouse in subset_dict:
        tot_ses   = subset_dict[mouse]['conc']['b_meta']['nsessions'] # total number of sessions for a given mouse
        ses_count = list(range(1,tot_ses+1))
        sessions  = list(set(ses_count).intersection(ses_n))
        if sessions: # if not an empty list - ie if the mouse has done sessions in the range relevant, as given by ses_n
            print(mouse)
            for ses in sessions:
                lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(subset_dict,mouse,ses)  
                p_coefs,r2,perf_drop = p_linr_regr_all(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,subset_dict,lick_df,photo_df,vid_df,bmeta,pmeta,beh_predictors,photo_predictors,ses,tw_start,tw_length,rew_only,ridge_regression=True,scale=True,tobepred=tobepred)
                if ses == 1:
                    nac_pcoefs1.append(p_coefs), r2s1.append(r2),perf_drops1.append(perf_drop)
                elif ses == 2:
                    nac_pcoefs2.append(p_coefs), r2s2.append(r2),perf_drops2.append(perf_drop)
    
    nac_pcoefs = [nac_pcoefs1,nac_pcoefs2]
    # t test with benjamini hochberg correction
    if tobepred == 'full':
        adj_ps1, adj_ps2 = [],[]
        for ses in ses_n:
            by_pred        = [pd.concat([nac_pcoefs[ses-1][i].iloc[:,j] for i in range(len(nac_pcoefs[ses-1]))],axis=1) for j in range(len(nac_pcoefs[ses-1][0].columns))]
            for pred in by_pred:
                ppvals                     = [stats.ttest_1samp(np.array(pred.iloc[i]),0)[1] for i in range(len(pred))]
                bool_p, adjusted_pps, _, _ = multipletests(ppvals, method='fdr_bh') # p_adjusted is boolean 
                if ses == 1:
                    adj_ps1.append(adjusted_pps)
                elif ses == 2:
                    adj_ps2.append(adjusted_pps)

        
    if plot:
        if tobepred == 'full':
            plot_coefs_window_ses(nac_pcoefs1,nac_pcoefs2,adj_ps1,adj_ps2,0.01,predic_labels,tw_start,tw_length)
            
        elif tobepred == 'peak':
            plot_perf_drop(perf_drops1, perf_drops2)
    


def p_linr_regr_all(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,data_dict,lick_df,photo_df,vid_df,bmeta,pmeta,beh_predictors,photo_predictors,ses_n,tw_start=-2,tw_length=5,rew_only=True,ridge_regression=True,scale=False,tobepred='peak'):
    """
    Run a linear or ridge regression predicting photometry responses from behavioral predictors.
    
    Builds a predictor matrix for one session, aligns predictors to visit-related
    photometry windows, and predicts either the full response trace or a summary
    value of the response.
    
    Parameters
    ----------
    travels_filepath : str
        Path to saved travel-related data.
    div_from_rand_filepath : str
        Path to saved divergence-from-random data.
    mat_GLM_path : str
        Path to GLM/AQUA model data.
    GLM_prob_filepath : str
        Path to GLM probability data.
    data_dict : dict
        Behavioral/session data.
    lick_df : pandas.DataFrame
        Lick/visit data.
    photo_df : pandas.DataFrame
        Photometry data containing a 'signal' column.
    vid_df : pandas.DataFrame
        Video-derived data.
    bmeta : dict-like
        Behavioral metadata.
    pmeta : dict-like
        Photometry metadata.
    beh_predictors : list
        Behavioral predictors to keep before final predictor selection.
    photo_predictors : list
        Final predictor columns used in the regression.
    ses_n : int
        Session number.
    tw_start : float, default=-2
        Start of photometry window relative to visit.
    tw_length : float, default=5
        Length of photometry window.
    rew_only : bool, default=True
        Whether to include only rewarded visits.
    ridge_regression : bool, default=True
        Whether to use ridge regression with cross-validated alpha.
    scale : bool, default=False
        Whether to z-score predictors before ridge regression.
    tobepred : {'full', 'peak', 'mean', 'smooth'}, default='peak'
        Photometry target to predict.
    
    Returns
    -------
    tuple
        Regression coefficients, model R², and performance drops from removing
        each predictor.
    """
    pred_df         = rvp.gen_predictors(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,data_dict,lick_df,vid_df,bmeta,[ses_n],data_type='real',task='conc',import_travels=False,no_travel=True)
    obs             = rvp.gen_observed(lick_df,bmeta)
    
    # add extra predictors into pred_df
    pred_df['photo'] = obs['photo_i']
    ses_photo        = photo_df['signal'].to_numpy()
    
    if rew_only:
        rew_val = pred_df['rewarded'].max()
        pred_df = pred_df[pred_df['rewarded']==rew_val].reset_index(drop=True)
    else:
        # if not reward only, remove alpha/AQUA_rpe/Q_rpe
        if np.any(np.isin(['alpha','Q_rpe','AQUA_rpe'],pred_df.columns)):
            pred_df = pred_df[[col for col in pred_df.columns if col not in ['alpha','Q_rpe','AQUA_rpe']]]

    # remove nan rows
    pred_df = pred_df.dropna(axis=1,how='all') # drop columns that only have nan first 
    pred_df = pred_df.dropna()

    # get photometry data in window 
    wini,_  = sf.extract_window(pred_df['photo'],ses_photo,pmeta['sampling_rate'][0],tw_start=tw_start,tw_length=tw_length)
    sig     = np.vstack([ses_photo[wini[x][0]:wini[x][1]] for x in range(len(wini))])
    
    # predicting single value for reward response - peak response - rather than response over time
    peak_signal = [np.max(sig[i]) for i in range(len(sig))]
    mean_signal = [np.mean(sig[i]) for i in range(len(sig))]

    # remove photometry indices and other predictors not wanting to be included in the model
    pred_df = pred_df.drop(['photo'],axis=1)
    pred_df = pred_df[beh_predictors]
    pred_df = pred_df[photo_predictors].reset_index(drop=True)
    
    # assigning variables: x are predictors, y is target - to be predicted
    X = pred_df 
    if tobepred == 'full':
        y = sig
    elif tobepred == 'peak':
        y = peak_signal
    elif tobepred == 'mean':
        y = mean_signal
    elif tobepred == 'smooth':
        y = list(pd.DataFrame(peak_signal).rolling(window=1,min_periods=1).mean().iloc[:,0])

    if ridge_regression:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        
        if scale:
            scaler  = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test  = scaler.transform(X_test)

        alphas = np.logspace(-4, 4, 100)  # Example range
        ridge_cv = RidgeCV(alphas=alphas, cv=5)
        ridge_cv.fit(X_train, y_train)
        
        # Get the best alpha from cross-validation
        best_alpha = ridge_cv.alpha_

        # Fit the Ridge model with the best alpha
        ridge = Ridge(alpha=best_alpha)
        ridge.fit(X_train, y_train)

        # calculate r2
        r2 = ridge.score(X_test,y_test)

        pcoefs = pd.DataFrame(ridge.coef_)
        if tobepred=='full':
            pcoefs.columns = X.columns

        from sklearn.model_selection import cross_val_score
        
        # Evaluate model performance with all predictors
        scores_all = cross_val_score(ridge, X_train, y_train, cv=5, scoring='r2')

        # Evaluate model performance with one predictor removed
        perf_drops = []
        for i in range(len(X_train[0])):
            X_train_scaled_reduced = np.delete(X_train, i, axis=1)
            scores_reduced         = cross_val_score(ridge, X_train_scaled_reduced, y_train, cv=5, scoring='r2')
            performance_drop = np.mean(scores_all) - np.mean(scores_reduced)
            perf_drops.append(performance_drop)

    else:
        # run standard linear regression over all data 
        regr           = linear_model.LinearRegression()
        run_regr       = regr.fit(X,y)
        pcoefs         = pd.DataFrame(run_regr.coef_)
        if tobepred=='full':
            pcoefs.columns = photo_predictors
        r2             = run_regr.score(X,y)
    
        
    return pcoefs,r2,perf_drops

### PLOTTING

def plot_coefs_window_ses(pcoefs1,pcoefs2,adj_ps1,adj_ps2,pval_thresh,predic_labels,tw_start,tw_length):
    """
    Plot time-resolved regression coefficients for selected predictors.
    
    For each predictor, plots the mean coefficient trace across animals for two
    sessions, with SEM shading and significance markers based on adjusted
    p-values.
    
    Parameters
    ----------
    pcoefs1 : list of pandas.DataFrame
        Time-resolved coefficient DataFrames for session 1.
    pcoefs2 : list of pandas.DataFrame
        Time-resolved coefficient DataFrames for session 2.
    adj_ps1 : list
        Adjusted p-values for session 1, ordered to match coefficient columns.
    adj_ps2 : list
        Adjusted p-values for session 2, ordered to match coefficient columns.
    pval_thresh : float
        Significance threshold for adjusted p-values.
    predic_labels : list
        Predictor names to plot.
    tw_start : float
        Window start time relative to event.
    tw_length : float
        Window length in seconds.
    """
    # subset each dataframe by predictor names
    pcoefs1_sub = [df[predic_labels] for df in pcoefs1]
    pcoefs2_sub = [df[predic_labels] for df in pcoefs2]
    
    # get indices of those predictors (based on column order)
    predic_idx = [pcoefs1[0].columns.get_loc(label) for label in predic_labels]
    
    # use indices to subset the corresponding entries in adj_ps1
    adj_ps1_sub = [adj_ps1[i] for i in predic_idx]
    adj_ps2_sub = [adj_ps2[i] for i in predic_idx]
        
    cmap    = load_cmap("Starfish")
    mcolors = [cmap(i) for i in [0,0.25,0.5,0.75,1.0]]
    
    # turning pvals into bool 
    bool_ps1 = [adj_ps1_sub[i]<pval_thresh for i in range(len(adj_ps1_sub))]
    bool_ps2 = [adj_ps2_sub[i]<pval_thresh for i in range(len(adj_ps2_sub))]
    
    pred_means1 = pd.concat([pcoefs1_sub[i] for i in range(0,len(pcoefs1_sub))]).groupby(level=0).mean()
    pred_sems1  = pd.concat([pcoefs1_sub[i] for i in range(0,len(pcoefs1_sub))]).groupby(level=0).sem()
    pos_sems1   = pred_means1+pred_sems1
    neg_sems1   = pred_means1-pred_sems1
    
    pred_means2 = pd.concat([pcoefs2_sub[i] for i in range(0,len(pcoefs2_sub))]).groupby(level=0).mean()
    pred_sems2  = pd.concat([pcoefs2_sub[i] for i in range(0,len(pcoefs2_sub))]).groupby(level=0).sem()
    pos_sems2   = pred_means2+pred_sems2
    neg_sems2   = pred_means2-pred_sems2

    for i in range(len(pred_means1.iloc[0,:])):
        fig,ax = plt.subplots(figsize=(2.3,2))
    
        ax.set_xlim([0,5])
        x       = range(0,len(pred_means1.iloc[:,i]))
        x_sec   = [sample/130 for sample in x]
        x_old   = np.array([0,-tw_start,tw_length])
        x_ticks = np.array([tw_start,0,tw_length+tw_start])
        ax.set_xticks(x_old,x_ticks)

        [ax.spines[loc].set_visible(False) for loc in ['top','right','left']]
        ax.plot(x_sec,pred_means1.iloc[:,i],color=mcolors[0])
        ax.plot(x_sec,pred_means2.iloc[:,i],color=mcolors[3])
        ax.fill_between(x_sec,pos_sems1.iloc[:,i],neg_sems1.iloc[:,i],alpha=0.2,color=mcolors[0])
        ax.fill_between(x_sec,pos_sems2.iloc[:,i],neg_sems2.iloc[:,i],alpha=0.2,color=mcolors[3])
        
        ax.axhline(y=0,xmin=0.1,color='gainsboro',linestyle=(0,(1,1)),linewidth=2,zorder=-1)
        ax.axvline(x=abs(tw_start),color='gainsboro',linewidth=2,zorder=-1)
        
        # p val line. true to false = -1. False to true = 1. 
        if np.any(bool_ps1[i]):
            ax.axhline(y=0.8,xmin=np.argmax(bool_ps1[i])/len(bool_ps1[i]),xmax=np.where(bool_ps1[i])[0][-1]/len(bool_ps1[i]),linewidth=3,color=mcolors[0])
        if np.any(bool_ps2[i]):
            ax.axhline(y=0.84,xmin=np.argmax(bool_ps2[i])/len(bool_ps2[i]),xmax=np.where(bool_ps2[i])[0][-1]/len(bool_ps2[i]),linewidth=3,color=mcolors[3])
        
        ax.set_ylim([-0.2,1])
        ax.set_yticks([])
        ax.set_yticklabels([])
        
        ax.yaxis.set_tick_params(width=1.5)
        ax.xaxis.set_tick_params(width=1.5)
        ax.tick_params(axis='x',labelsize=18,direction='in')
        ax.tick_params(axis='y',labelsize=18,direction='out')

        ax.spines[['bottom']].set_linewidth(1.5)

        fig.tight_layout()


def plot_perf_drop(perf_drops1, perf_drops2):
    """
    Plot predictor performance drops for two sessions.
    
    Creates interleaved boxplots showing the cross-validated reduction in R² when
    each predictor is removed, with one box per predictor per session.
    
    Parameters
    ----------
    perf_drops1 : list
        Performance drops for session 1, with one list of predictor drops per
        animal.
    perf_drops2 : list
        Performance drops for session 2, with one list of predictor drops per
        animal.
    """

    hfont = {'fontname': 'Arial'}

    df1 = pd.DataFrame(perf_drops1)
    df2 = pd.DataFrame(perf_drops2)

    fig, ax = plt.subplots(figsize=(6, 3))

    cmap    = load_cmap("Starfish")
    mcolors = [cmap(i) for i in [0, 0.25, 0.5, 0.75, 1.0]]

    x1 = np.arange(0, 16, 2)
    x2 = x1 + 1

    for i in range(df1.shape[1]):

        # session 1 (e.g. mcolors[0])
        sns.boxplot(
            y=df1.iloc[:, i], x=[x1[i]] * len(df1),
            ax=ax, color=mcolors[0], fill=False,
            width=0.4, showfliers=False, showcaps=False
        )
        sns.stripplot(
            y=df1.iloc[:, i], x=[x1[i]] * len(df1),
            ax=ax, edgecolor='gainsboro', facecolor='white',
            linewidth=2, zorder=-1
        )

        # session 2 (e.g. mcolors[3])
        sns.boxplot(
            y=df2.iloc[:, i], x=[x2[i]] * len(df2),
            ax=ax, color=mcolors[3], fill=False,
            width=0.4, showfliers=False, showcaps=False
        )
        sns.stripplot(
            y=df2.iloc[:, i], x=[x2[i]] * len(df2),
            ax=ax, edgecolor='gainsboro', facecolor='white',
            linewidth=2, zorder=-1
        )

    ax.set_ylabel('C.V. reduction in r2', **hfont, fontsize=18)
    ax.set_ylim([-0.02, 0.15])
    ax.set_yticks([0, 0.15])
    ax.set_yticklabels([0, 0.15], fontsize=18)
    ax.set_xlim([-1, 16])

    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)

    ax.axhline(y=0, xmin=0.1, color='gainsboro',
               linestyle=(0,(1,1)), linewidth=2, zorder=-1)

    fig.tight_layout()

