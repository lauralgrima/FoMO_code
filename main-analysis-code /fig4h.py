import import_mat as im
import matplotlib.pyplot as plt
import seaborn as sns


def DA_as_alpha(comp_filepath):
    '''
    Plots for AQUA as DA - fig 7h and 7i
    '''
    r2,_,delta_rew,_ = im.import_modelcomp(comp_filepath)

    for i,metric in enumerate([r2,delta_rew]):
        hfont    = {'fontname':'Arial'} 
        fig,ax   = plt.subplots(figsize=(3,3))
            
        sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    
        sns.stripplot(y=metric['AQUA opt com'],x=0,edgecolor='silver',facecolor='white',linewidth=2,zorder=-1)
        sns.stripplot(y=metric['AQUA DA=alpha'],x=0.2,edgecolor='silver',facecolor='white',linewidth=2,zorder=-1)
        
        sns.boxplot(y=metric['AQUA opt com'],x=0.2,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
        sns.boxplot(y=metric['AQUA DA=alpha'],x=0,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
        
        [sns.lineplot(x=[0,1],y=[metric['AQUA opt com'][i],metric['AQUA DA=alpha'][i]],ax=ax,color='silver',linewidth=2,zorder=-2) for i in range(0,len(metric['AQUA opt com']))]
    
        if i == 0:
            ax.set_ylabel(r'trans. matrix similarity (r$^2$)',**hfont,fontsize=18)
            
            ax.set_yticks([0,0.5,1])
            ax.set_yticklabels([0,0.5,1],fontsize=18)
            ax.set_ylim([0,1])
        
        else:
            ax.set_ylabel(r'$\Delta$ predicted rews. col.',**hfont,fontsize=18)
            
            ax.set_yticks([-200,0,100])
            ax.set_yticklabels([-200,0,100],fontsize=18)
            ax.set_ylim([-200,100])
        
        ax.set_xlim([-1,5])
        
        ax.spines[['bottom','top','right']].set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.yaxis.set_tick_params(width=1.5)
    
        fig.tight_layout()
