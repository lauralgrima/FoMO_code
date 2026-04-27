import import_mat as im
import support_funcs as sf
import seaborn as sns
import matplotlib.pyplot as plt


def MULTImat_similarity(data_dict,mat_GLM_path,ses_n=1,plot=True):
    '''
    Extract AQUA model r² values across mice for a single session.

    Parameters:
        data_dict: dict
            Full behavioral data dictionary.

        mat_GLM_path: str
            Path to AQUA/GLM .mat data.

        ses_n: int
            Session number to extract. Uses 1-indexing.

        plot: bool
            If True, plot the distribution of r² values.

    Returns:
        r2s: list
            AQUA model r² values for mice with data for session ses_n.
    '''
    # extracting just interval task, non opto mice, config 1
    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    subset_dict     = {**subset_dict_NAc,**subset_dict_DMS}
    
    # importing AQUA model data 
    _,AQUA_beh_NAc = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_beh = {**AQUA_beh_DMS,**AQUA_beh_NAc}
    
    r2s = [AQUA_beh[mouse]['r2'][ses_n-1] for mouse in list(subset_dict.keys()) if len(AQUA_beh[mouse]['r2'])>= ses_n]
    
    if plot:
        plot_mat_similarity(r2s)
    
    return r2s


def plot_mat_similarity(r2s):
    '''
    Plot distribution of transition matrix similarity (r values).
    '''

    hfont    = {'fontname':'Arial'} 
    fig,ax   = plt.subplots(figsize=(3,3))

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    sns.boxplot(y=r2s,x=0,ax=ax,color='black',fill=False,width=0.5,showfliers=False,showcaps=False)
    sns.stripplot(y=r2s,x=0,edgecolor='silver',facecolor='white',linewidth=2,zorder=-1)
    
    ax.set_ylabel(r'trans. matrix similarity ($r$)',**hfont,fontsize=18)
    
    ax.set_yticks([0,0.5,1])
    ax.set_yticklabels([0,0.5,1],fontsize=18)
    ax.set_xlim([-1,5])
    ax.set_ylim([0,1])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)

    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()
