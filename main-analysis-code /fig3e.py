import fig2b
import support_funcs as sf
import import_mat as im
import matplotlib.pyplot as plt
import seaborn as sns


def MULTImatch_AQUA(data_dict,mat_GLM_path,ses_n=1,plot=True):
    '''
    Compute matching behavior (reward vs. choice ratios) across mice using AQUA model data.

    For each mouse in the interval task (config 1, NAc and DMS), this function extracts
    AQUA-predicted visits and computes reward and choice (visit) ratios using
    `matching_calc`. Results are aggregated across mice.

    Parameters
    ----------
    data_dict : dict
        Full behavioral data dictionary.

    mat_GLM_path : str
        Path to AQUA/GLM .mat data.

    ses_n : int, optional
        Session number to analyze (1-indexed).

    plot : bool, optional
        If True, plot reward vs. choice ratios using `plot_match`.

    Returns
    -------
    amatch_rew : list of array-like
        Reward ratios for each mouse.

    amatch_vis : list of array-like
        Choice (visit) ratios for each mouse.

    Notes
    -----
    - Only mice with at least `ses_n` sessions are included.
    - Matching is computed in AQUA space (model-predicted behavior).
    - Port visits are converted to rank space prior to analysis.
    '''
    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    subset_dict     = {**subset_dict_NAc,**subset_dict_DMS}
    
    _,AQUA_beh_NAc  = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS  = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_beh        = {**AQUA_beh_DMS,**AQUA_beh_NAc}
    
    amatch_rew,amatch_vis,aslopes = [],[],[]
    for mouse in subset_dict:
      print(mouse)
      tot_ses = data_dict[mouse]['conc']['b_meta']['nsessions']
      if ses_n <= tot_ses:
          lick_df,_,_,bmeta,_   = sf.extract_data(subset_dict,mouse,ses_n)
          AQUA_visits           = im.get_visits_AQUA(AQUA_beh[mouse]['AQUA_vis_invis'],lick_df,ses_n,bmeta)['port']
          AQUA_visits           = AQUA_visits.map(sf.port_2_rank(bmeta,ses_n))
          
          match_rew,match_vis,fit_parameters = fig2b.matching_calc(lick_df,bmeta,ses_n,plot_ind=False,data='AQUA')
          amatch_rew.append(match_rew),amatch_vis.append(match_vis)
          
          # sensitivity plot
          if len(fit_parameters) > 0:
              aslopes.append(fit_parameters[0][0])
    
    if plot: 
        fig2b.plot_match(amatch_rew,amatch_vis,data='AQUA')
        plot_AQUAsensi(aslopes)


### PLOTTING 

def plot_AQUAsensi(aslopes):
    '''
    Plot distribution of transition matrix similarity (r values).
    '''

    hfont    = {'fontname':'Arial'} 
    fig,ax   = plt.subplots(figsize=(3,3))

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    sns.boxplot(y=aslopes,x=0,ax=ax,color='black',fill=False,width=0.5,showfliers=False,showcaps=False)
    sns.stripplot(y=aslopes,x=0,edgecolor='silver',facecolor='white',linewidth=2,zorder=-1)
    
    ax.set_ylabel('sensitivity',**hfont,fontsize=18)
    
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
    ax.set_yticklabels([0,0.2,0.4,0.6,0.8,1],fontsize=18)
    ax.set_xlim([-1,5])
    ax.set_ylim([0,1])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)

    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()

