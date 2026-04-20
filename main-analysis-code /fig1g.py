import pandas as pd
import numpy as np
import support_funcs as sf
import strategy_benchmarks as sb
import matplotlib.pyplot as plt
        

def eg_rate_of_return(data_dict,mouse='6PG6',ses_n=1,win_size=2500,nstrategy_samples=100,plot=True):
    """
    Compare reward return over time for mouse and simulated benchmarks.

    Computes sliding-window reward rate (rewards/check) for the animal,
    random, ideal observer, and WSLS strategies for one mouse and session.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary of loaded data.
    mouse : str, default='6PG6'
        Subject ID to analyze.
    ses_n : int, default=1
        Session number.
    win_size : int, default=2500
        Sliding window size, in 1 s bins.
    nstrategy_samples : int, default=100
        Number of strategy simulations to average over.
    plot : bool, default=True
        If True, plot the sliding-window comparison.

    Returns
    -------
    None
    """
    lick_df  = data_dict[mouse]['conc']['mlick_df']
    bmeta    = data_dict[mouse]['conc']['b_meta']
    
    # mouse data 
    visit_times = np.array(lick_df.loc[lick_df['unique_visit']==1]['event_time'])
    rew_times   = np.array(lick_df.loc[(lick_df['rewarded']==1)&(lick_df['unique_visit']==1)]['event_time'])

    longform_visits  = sf.longform(visit_times) # in 1s bins
    longform_rewards = sf.longform(rew_times)
    
    av_checks_win  = pd.DataFrame(longform_visits).rolling(window=win_size,min_periods=1).mean()
    av_rews_win    = pd.DataFrame(longform_rewards).rolling(window=win_size,min_periods=1).mean()
    av_return_win  = av_rews_win/av_checks_win

    # simulated strategies
    rand_income_multisampled = []
    io_income_multisampled   = []
    wsls_income_multisampled = []
    for n in range(nstrategy_samples+1):
        
        # random strategy
        print(n)
        rand_strategy_df = sb.simulate_strategies(visit_times,bmeta,ses_n,strategy='random')
        rand_rew_times   = rand_strategy_df.loc[rand_strategy_df['outcomes']==1]['check_times']
        longform_rand    = sf.longform(rand_rew_times)
        
        rand_rews_win    = pd.DataFrame(longform_rand).rolling(window=win_size,min_periods=1).mean()
        rand_income_win  = rand_rews_win/av_checks_win
        
        rand_income_multisampled.append(rand_income_win)
    
        # ideal observer 
        io_strategy_df = sb.simulate_strategies(visit_times,bmeta,ses_n,strategy='io')
        io_rew_times   = io_strategy_df.loc[io_strategy_df['outcomes']==1]['check_times']
        longform_io    = sf.longform(io_rew_times)
        
        io_rews_win    = pd.DataFrame(longform_io).rolling(window=win_size,min_periods=1).mean()
        io_income_win  = io_rews_win/av_checks_win
        
        io_income_multisampled.append(io_income_win)
        
        # WSLS
        wsls_strategy_df = sb.simulate_strategies(visit_times,bmeta,ses_n,strategy='winstay')
        wsls_rew_times   = wsls_strategy_df.loc[wsls_strategy_df['outcomes']==1]['check_times']
        longform_wsls    = sf.longform(wsls_rew_times)
        
        wsls_rews_win    = pd.DataFrame(longform_wsls).rolling(window=win_size,min_periods=1).mean()
        wsls_income_win  = wsls_rews_win/av_checks_win
        
        wsls_income_multisampled.append(wsls_income_win)

    io_income_mean   = pd.concat(io_income_multisampled,axis=1).mean(axis=1)
    io_income_sem    = pd.concat(io_income_multisampled,axis=1).sem(axis=1)

    # make relative to IO
    prop_return = (av_return_win[0]/io_income_mean)
    rand_prop   = [(rand_income_multisampled[i][0]/io_income_mean) for i in range(len(rand_income_multisampled))]
    wsls_prop   = [(wsls_income_multisampled[i][0]/io_income_mean) for i in range(len(wsls_income_multisampled))]
    
    rand_prop_mean = pd.DataFrame(rand_prop).transpose().mean(axis=1)
    rand_prop_sem  = pd.DataFrame(rand_prop).transpose().sem(axis=1)
    wsls_prop_mean = pd.DataFrame(wsls_prop).transpose().mean(axis=1)
    wsls_prop_sem  = pd.DataFrame(wsls_prop).transpose().sem(axis=1)
    

    if plot:
        plotting_eg_rate_of_return(prop_return,rand_prop_mean,rand_prop_sem,io_income_mean,io_income_sem,wsls_prop_mean,wsls_prop_sem)


### PLOTTING

def plotting_eg_rate_of_return(av_income_win,rand_income_mean,rand_income_sem,io_income_mean,io_income_sem,wsls_income_mean,wsls_income_sem):
    '''
    av_income_win and io_income_win should be one dimensional numpy arrays. 
    rand_income_multisampled should be a list of length nstrategy_samples, of one dimensional numpy arrays.
    set plot_prop to True to plot proportion of IO values rather than rate of return per se
    '''

    # set up
    hfont = {'fontname':'Arial'} 
    fig, ax = plt.subplots(figsize=(3.5,3.5))
    
    # plot
    ax.plot(range(len(av_income_win)),av_income_win,color='dimgrey')
    ax.plot(range(len(rand_income_mean)),rand_income_mean,color='orange',zorder=-1)
    ax.plot(range(len(wsls_income_mean)),wsls_income_mean,color='peru',zorder=-1)

    ax.set_ylim([0.5,1])
    ax.set_yticks([0.5,1.0])
    ax.set_xlim([0,10800])
    ax.set_xticks([0,60*60,120*60,180*60])
    ax.set_xticklabels([0,60,120,180],**hfont,fontsize=18)
    
    ax.spines[['right','top']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both', which='major', labelsize=18)
    
    fig.tight_layout()
    plt.show()

