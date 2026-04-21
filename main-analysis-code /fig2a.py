import numpy as np
import support_funcs as sf
import matplotlib.pyplot as plt

def cumu_visits(data_dict,mouse='6PG6',ses_n=1,sort=True,plot=True):
    '''
    Exract times and numbers of cumulative visits for each port for a session.
    '''
    lick_df  = data_dict[mouse]['conc']['mlick_df']
    bmeta    = data_dict[mouse]['conc']['b_meta']
    
    port_times = [lick_df[(lick_df['unique_visit']==1)&(lick_df['port']==port)]['event_time'] for port in range(1,7)]
    cumu_vis   = [np.arange(1,len(port_times[i])+1) for i in range(0,6)]
      
    if sort: # sort by rank
        port_times = sf.sort_data_by_interval(port_times,bmeta['intervals'][ses_n-1])
        cumu_vis   = sf.sort_data_by_interval(cumu_vis,bmeta['intervals'][ses_n-1])
          
    if plot: # plot for a single session
        plot_cumu_visits(port_times,cumu_vis)
        
          
def plot_cumu_visits(port_times,cumu_vis):
    '''
    Plot cumulative visit curves split by port for a single session.
    '''
    
    fig,ax = plt.subplots(figsize=(3.5,3.5))
    hfont  = {'fontname':'Arial'}
    color  = plt.cm.plasma(np.linspace(0,0.9,6))
    
    for port in range(0,6):
        ax.plot(port_times[port]/60,cumu_vis[port],color=color[port],label=port+1)
        
    ax.set_xlabel('time (min.)',**hfont,fontsize=18)
    ax.set_ylabel('cumulative visits',**hfont,fontsize=18)
    ax.spines[['right','top']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.set_xticks([0,60,120,180])
    ax.set_yticks([0,350,700])
    ax.set_xticklabels([0,60,120,180],**hfont,fontsize=18)
    ax.set_yticklabels([0,350,700],**hfont,fontsize=18)
    ax.set_xlim([0,180])
    ax.set_ylim([0,700])
    
    fig.tight_layout()
