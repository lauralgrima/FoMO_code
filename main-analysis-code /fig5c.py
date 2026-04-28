import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd


def opto_calibration(opto_stim_folder,tw_start=-2,tw_length=5,plot=True):
    """
    Plot mean photometry responses for optogenetic stimulation and reward.
    
    Loads precomputed per-animal response traces from CSV files, computes the
    across-animal mean and SEM for stimulation and reward conditions, and
    optionally plots the resulting traces.
    
    Parameters
    ----------
    opto_stim_folder : str
        Path to folder containing 'opto_stim_df.csv' and 'opto_rew_df.csv'.
    tw_start : float, default=-2
        Window start time relative to event.
    tw_length : float, default=5
        Window length in seconds.
    plot : bool, default=True
        Whether to generate plots of the mean ± SEM traces.
    
    Returns
    -------
    None
    """

    # import stim and rew dataframes
    anm_stim_mean = pd.read_csv(opto_stim_folder+'opto_stim_df.csv') # rows are mice, columns are timepoints
    anm_rew_mean  = pd.read_csv(opto_stim_folder+'opto_rew_df.csv') # rows are mice, columns are timepoints
    
    anm_stim_mean = anm_stim_mean.values.tolist()
    anm_rew_mean  = anm_rew_mean.values.tolist()

    all_stim_mean = np.mean(anm_stim_mean,axis=0)
    all_rew_mean = np.mean(anm_rew_mean,axis=0)
    
    all_stim_sem               = stats.sem(anm_stim_mean,axis=0,nan_policy='omit')
    all_rew_sem                = stats.sem(anm_rew_mean,axis=0,nan_policy='omit')
    all_stim_pos               = all_stim_mean+all_stim_sem
    all_stim_neg               = all_stim_mean-all_stim_sem
    all_rew_pos                = all_rew_mean+all_rew_sem
    all_rew_neg                = all_rew_mean-all_rew_sem
    
    if plot:
        plot_trace(tw_start,tw_length,all_stim_mean, all_stim_pos, all_stim_neg, 'deepskyblue', 'stim')
        plot_trace(tw_start,tw_length,all_rew_mean, all_rew_pos, all_rew_neg, 'navy', 'reward')


def plot_trace(tw_start,tw_length,mean,pos,neg,color,label):
    """
    Plot a mean time-series trace with shaded error.
    
    Plots a time-aligned signal (e.g., photometry) with upper and lower bounds
    (e.g., mean ± SEM), including basic formatting and reference lines.
    
    Parameters
    ----------
    tw_start : float
        Window start time relative to event.
    tw_length : float
        Window length in seconds.
    mean : array-like
        Mean signal trace.
    pos : array-like
        Upper bound of the signal (e.g., mean + SEM).
    neg : array-like
        Lower bound of the signal (e.g., mean - SEM).
    color : str or tuple
        Line and fill color.
    label : str
        Label for the trace.
    
    Returns
    -------
    None
    """
    
    hfont = {'fontname': 'Helvetica'} 
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    x = list(range(len(mean)))
    x_sec = [sample / 130 for sample in x]
    x_old = [0, 0 - tw_start, tw_length]
    x_ticks = [tw_start, 0, tw_length + tw_start]

    ax.plot(x_sec, mean, color=color, label=label)
    ax.fill_between(x_sec, pos, neg, alpha=0.3, color=color)

    ax.axvline(x=abs(tw_start), color='gainsboro', linewidth=2, zorder=-1)
    ax.axhline(y=0, color='gainsboro', linestyle=(0, (1, 1)), linewidth=2, zorder=-1)

    ax.set_xticks(x_old)
    ax.set_xticklabels(x_ticks)
    ax.set_ylim([-0.5, 3])
    ax.set_yticks([0, 1])
    ax.set_xlim([0, 5])

    ax.spines[['right', 'top', 'left']].set_visible(False)
    ax.spines[['bottom']].set_position(('outward', 10))
    ax.tick_params(axis='both', which='both', labelsize=18, direction='in')
    ax.yaxis.set_tick_params(width=1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.spines[['left', 'bottom']].set_linewidth(1.5)

    ax.set_xlabel('time (s)', **hfont, fontsize=18)

    fig.tight_layout()
    plt.show()

    




