import import_mat as im
import matplotlib.pyplot as plt
import seaborn as sns

def MULTImodel_comp(comp_filepath,plot_type='comp_matrix'):
    
    r2,inc,delta_rew,_    =  im.import_modelcomp(comp_filepath)
    
    if plot_type == 'comp_matrix':
        data = r2
    elif plot_type == 'rew_col':
        data = delta_rew
    
    # select appropriate columns 
    restr_data = data.iloc[:,[0,1,5,6,7,8]]
    
    hfont   = {'fontname':'Arial'} 
    fig,ax  = plt.subplots(figsize=(4,3))

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    x = [0,1,4,5,6,7]
    for i in range(len(restr_data.iloc[0,:])):
        sns.boxplot(y=restr_data.iloc[:,i],x=x[i],ax=ax,color='black',fill=False,width=0.3,showfliers=False,showcaps=False,)
        sns.stripplot(y=restr_data.iloc[:,i],x=x[i],edgecolor='gainsboro',facecolor='white',linewidth=2,zorder=-1)

    ax.set_ylabel(r'$\Delta$ rewards collected',**hfont,fontsize=18)
    
    if plot_type == 'rew_col':
        ax.set_ylabel(r'$\Delta$ rewards collected',**hfont,fontsize=18)
        ax.set_ylim([-200,100])
        ax.set_yticks([-200,-100,0,100])
        ax.set_yticklabels([-200,-100,0,100],fontsize=18)
    else:
        ax.set_ylabel('trans. matrix similarity (r)',**hfont,fontsize=18)
        ax.set_ylim([-0.2,1])
        ax.set_yticks([0,0.5,1])
        ax.set_yticklabels([0,0.5,1],fontsize=18)
       
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
 
    ax.axhline(y=0,xmin=0.1,color='gainsboro',linestyle=(0,(1,1)),linewidth=2,zorder=-1)
       
    fig.tight_layout()  
