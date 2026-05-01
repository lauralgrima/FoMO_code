import load_all_data as lad
import fig1c, fig1e, fig1f, fig1g, fig1h
import fig2a, fig2b, fig2c, fig2d, fig2e, fig2g
import fig3b, fig3c, fig3d, fig3e, fig3f, fig3g, fig3hi
import fig4b, fig4c, fig4df, fig4e, fig4h
import fig5a, fig5c, fig5d
import fig6bc, fig6d, fig6e, fig6f, fig6g

# for information about which experiments mice did etc. see the 'data_summary' word doc uploaded in the github repo

subject_IDs = ['6PG5','6PG6','6PG8','6PG9','6PG10','6PG11','6PG12','6PG15','6PG24','6PG25','6PG27','6PG28','6PG30','6PG31','6PG33','6PG34','6PG35','6PO3','6PO4','6PO5',
               '6PO7','6PO8','6PO9','6PO10','dProb1','dProb3','gProb1','gProb2','ML11','ML12','ML13','ML14','ML15','ML17']
regions     = ['NAc','DMS']
tasks       = ['conc','prob'] # conc refers to the interval schedule task, prob refers to probabality schedule task 

data_dict   = lad.load_data(df_folder, subject_IDs, tasks)


### ----FIGURE 1---- ###

fig1c.eg_event_traces(data_dict, mouse='6PG6', cutoff=30, axes_off=True)
fig1e.MULTIdiverge_time(data_dict,ses_n=1,nstrategy_samples=100,slope_thresh=0.02,win_len=10,ind_plot=False,plot=True)
fig1f.MULTIperc_ideal(data_dict,ses_n=1,nstrategy_runs=100,plot=True)
fig1g.eg_rate_of_return(data_dict,mouse='6PG6',ses_n=1,win_size=2500,nstrategy_samples=100,plot=True)
fig1h.MULTIkl_diverge_ratio(data_dict,ses_n=1,final_policy='IO',consec_counts=100,nrand_runs=100,plot_ind=False,plot=True) # set ses_n to 2 to visualise session 2. set plot_ind to True for individual examples, as in left panel of 1h

### ----FIGURE 2---- ###

fig2a.cumu_visits(data_dict,mouse='6PG6',ses_n=1,sort=True,plot=True)
fig2b.eg_matching(data_dict,mouse='6PG6',ses_n=1)
fig2b.MULTImatching(data_dict,ses_n=1,plot=True)
fig2c.eg_learning_to_match(data_dict,mouse='6PG6',ses_n=1)
fig2c.MULTIlearning_to_match(data_dict,ses_n=1,win_len=30,plot=True)
fig2d.MULTIsensi_across_days(data_dict, ses_n=[1,2],plot=True)
fig2e.MULTIivis(data_dict,ses_n=1,spatial_config=1,plot=True)
fig2g.MULTIconditional_matching(data_dict,ses_n=1,plot=True)

### ----FIGURE 3---- ###

fig3b.eg_transition_matrix(data_dict,mat_GLM_path,mouse='6PG6',port_or_rank='rank',ses_n=1,plot=True)
fig3c.MULTImat_similarity(data_dict,mat_GLM_path,ses_n=1,plot=True)
fig3d.eg_cumu_visits_AQUA(data_dict,mat_GLM_path,mouse='6PG6',ses_n=1,plot=True)
fig3e.MULTImatch_AQUA(data_dict,mat_GLM_path,ses_n=1,plot=True)
fig3f.MULTI_kl_withinses(data_dict,mat_GLM_path,ses_n=1,win_len=50,max_vis=300,ses_fraction=4,portion='end',plot=True)
fig3g.MULTImultinomial_regr(data_dict,mat_GLM_path,ses_n=1,cv=5,pred_approach='train_test',plot=True)
fig3hi.MULTImodel_comp(comp_filepath,plot_type='comp_matrix')
fig3hi.MULTImodel_comp(comp_filepath,plot_type='rew_col')

### ----FIGURE 4---- ###

fig4b.MULTIvisit_DA(data_dict,ses_n=[1,2],tw_start=-2,tw_length=5,plot=True)
fig4c.MULTIvisit_DA_by_port(data_dict,ses_n=1,tw_start=-2,tw_length=5,plot=True)
fig4df.MULTIplinr_whole(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,data_dict,ses_n=[1,2],tw_start=-2,tw_length=5,plot=True,rew_only=False,tobepred='full')
fig4e.eg_ind_resp(data_dict,mouse='6PG6',tw_start=0,tw_length=2,ses_n=1,plot=True)
fig4df.MULTIplinr_whole(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,data_dict,photo_predictors,beh_predictors,ses_n=[1,2],tw_start=-2,tw_length=5,plot=True,rew_only=True,tobepred='full')
fig4df.MULTIplinr_whole(travels_filepath,div_from_rand_filepath,mat_GLM_path,GLM_prob_filepath,data_dict,photo_predictors,beh_predictors,ses_n=[1,2],tw_start=-2,tw_length=5,plot=True,rew_only=True,tobepred='peak')
fig4h.DA_as_alpha(comp_filepath)

### ----FIGURE 5---- ###

fig5a.opto_updaterule(opto_filepath,plot=True)
fig5c.opto_calibration(opto_filepath,tw_start=-2,tw_length=5,plot=True)
fig5d.MULTImatch_opto_4worst(data_dict,limit=90,ses_n=1,plot=True)

### ----FIGURE 6---- ###

fig6bc.MULTIses_pchoice(data_dict,mat_GLM_path,win_len=20,plot=True,plot_example=False)
fig6bc.MULTIses_pchoice(data_dict,mat_GLM_path,win_len=20,plot=False,plot_example=True)
fig6d.MULTImultises_kl(data_dict,mat_GLM_path,lim2=200,lim3=300,win_len=50,plot=True)
fig6e.MULTIses_rewpeak_DA(data_dict,tw_start=0,tw_length=2,n_quants=4)
fig6f.MULTIses3_DA(data_dict,tw_start=-2,tw_length=5,plot=True)
fig6g.MULTIalpha_KL(KL_path)
fig6g.MULTIDA_KL(KL_path,incl_ses1=False)