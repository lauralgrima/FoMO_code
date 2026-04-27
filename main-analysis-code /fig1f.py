import numpy as np
import support_funcs as sf
import strategy_benchmarks as sb
import matplotlib.pyplot as plt
import seaborn as sns


def MULTIperc_ideal(data_dict,ses_n=1,nstrategy_runs=100,plot=True):
    """
    Compare reward collection efficiency across mice.

    Computes the percentage of possible rewards collected (relative to an
    ideal observer) for random, WSLS, and animal behavior, for one session.
    Subsets to non-opto mice with the 'conc' task.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary of loaded data.
    ses_n : int
        Session number to analyze.
    nstrategy_runs : int, default=100
        Number of simulations for each strategy.
    plot : bool, default=True
        If True, plot summary results.

    Returns
    -------
    rand : list
        Percent reward (relative to ideal) for random strategy.
    WSLS : list
        Percent reward for win-stay-lose-switch strategy.
    anms : list
        Percent reward for animal behavior.
    """
    # extracting just interval task, non opto mice 

    subset_dict = sf.subset_mice(data_dict, task='conc', include_opto=False, config=None, region=None)
    
    rand,WSLS,anms = [],[],[]
    for mouse in subset_dict.keys:
        print(mouse)
        if subset_dict[mouse]['conc']['b_meta']['nsessions'] >= ses_n:
            lick_df,_,_,bmeta,_   = sf.extract_data(subset_dict,mouse,ses_n)  
            av_rand,av_WSLS,anm   = perc_ideal(lick_df,ses_n,bmeta,nstrategy_runs)
            rand.append(av_rand), WSLS.append(av_WSLS), anms.append(anm)

    if plot:
        plot_perc_ideal(rand,WSLS,anms)


def perc_ideal(lick_df,ses_n,bmeta,nstrategy_runs):
    """
    Compute reward collection efficiency relative to an ideal observer.

    Simulates ideal observer (IO), random, and WSLS strategies for a session,
    and compares their total rewards to the animal’s performance.

    Parameters
    ----------
    lick_df : pandas.DataFrame
        Lick/visit dataframe for one session.
    ses_n : int
        Session number.
    bmeta : dict
        Behavioral metadata used for simulation.
    nstrategy_runs : int
        Number of simulations for random and WSLS strategies.

    Returns
    -------
    av_rand : float
        Mean fraction of IO rewards collected by the random strategy.
    av_WSLS : float
        Mean fraction of IO rewards collected by the WSLS strategy.
    anm : float
        Fraction of IO rewards collected by the animal.
    """
    # mouse data 
    visit_times = np.array(lick_df.loc[lick_df['unique_visit']==1]['event_time'])
    tot_rews    = len(np.array(lick_df.loc[(lick_df['rewarded']==1)&(lick_df['unique_visit']==1)]['event_time']))
    
    # simulated strategies
    rand_percs,WSLS_percs = [],[]
    for n in range(nstrategy_runs+1):
        print(n)
        
        # OIO
        io_strategy_df = sb.simulate_strategies(visit_times,bmeta,ses_n,stagger=True,strategy='io')
        io_tot_rews    = io_strategy_df['outcomes'].sum()
        
        # random
        rand_strategy_df = sb.simulate_strategies(visit_times,bmeta,ses_n,stagger=True,strategy='random')
        rand_tot_rews    = rand_strategy_df['outcomes'].sum()
        rand_percs.append(rand_tot_rews/io_tot_rews)
        
        # WSLS
        WSLS_strategy_df = sb.simulate_strategies(visit_times,bmeta,ses_n,stagger=True,strategy='winstay')
        WSLS_tot_rews    = WSLS_strategy_df['outcomes'].sum()
        WSLS_percs.append(WSLS_tot_rews/io_tot_rews)
        
    av_rand = np.mean(rand_percs)
    av_WSLS = np.mean(WSLS_percs)
    anm     = tot_rews/io_tot_rews
    
    return av_rand,av_WSLS,anm


### PLOTTING

def plot_perc_ideal(rand,WSLS,anms):
    """
    Plot reward efficiency relative to ideal across strategies.

    Displays individual data points and boxplots for random, WSLS, and
    animal performance, with lines connecting matched subjects.

    Parameters
    ----------
    rand : list or array-like
        Percent of ideal rewards for random strategy.
    WSLS : list or array-like
        Percent of ideal rewards for WSLS strategy.
    anms : list or array-like
        Percent of ideal rewards for animals.

    Returns
    -------
    None
    """
    
    hfont = {'fontname':'Arial'} 
    fig,ax = plt.subplots(figsize=(3.5,3))
        
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    
    sns.stripplot(y=anms,x=0,edgecolor='gainsboro',facecolor='white',linewidth=2,zorder=-1)
    sns.stripplot(y=rand,x=0.2,edgecolor='gainsboro',facecolor='white',linewidth=2,zorder=-1)
    sns.stripplot(y=WSLS,x=0.4,edgecolor='gainsboro',facecolor='white',linewidth=2,zorder=-1)

    sns.boxplot(y=anms,x=0,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
    sns.boxplot(y=rand,x=0.2,ax=ax,color='orange',fill=False,width=0.4,showfliers=False,showcaps=False)
    sns.boxplot(y=WSLS,x=0.4,ax=ax,color='peru',fill=False,width=0.4,showfliers=False,showcaps=False)
    
    [sns.lineplot(x=[0,1,2],y=[anms[i],rand[i],WSLS[i]],ax=ax,color='gainsboro',linewidth=2,zorder=-2) for i in range(0,len(anms))]

    ax.set_ylabel('proportion of ideal',**hfont,fontsize=18)
    
    ax.set_yticks([0.4,0.6,0.8,1.0])
    ax.set_yticklabels([0.4,0.6,0.8,1.0],fontsize=18)
    ax.set_xlim([-1,5.5])
    ax.set_ylim([0.4,1.0])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()
    
