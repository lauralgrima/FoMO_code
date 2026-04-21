import support_funcs as sf
import matplotlib.pyplot as plt
import seaborn as sns
from fig2b import matching_calc


def MULTIsensi_across_days(data_dict, ses_n=[1, 2], plot=True):
    """
    Compare matching sensitivity across multiple sessions for each mouse.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary containing behavioral data for all subjects.
    ses_n : list of int, optional
        Session numbers to compare (default is [1, 2]).
    plot : bool, optional
        If True, plot session-to-session sensitivity values (default is True).

    Returns
    -------
    all_slopes : list of list
        Matching slopes for each included mouse across the requested sessions.
    config_type : list of int
        Configuration label for each included mouse.

    Notes
    -----
    Restricts analysis to subjects with `'conc'` data and excludes optogenetic
    mice (IDs starting with `'6PO'`). Only mice with all requested sessions are
    included. Configuration type is determined from the first session's interval
    structure in the raw metadata.
    """

    # extracting just interval task, non opto mice
    subset_dict = {
        subj: {'conc': data['conc']}
        for subj, data in data_dict.items()
        if 'conc' in data and not subj.startswith('6PO')
    }

    all_slopes = []
    config_type = []

    required_sessions = set(ses_n)

    for mouse in subset_dict.keys():
        print(mouse)

        nsessions = subset_dict[mouse]['conc']['b_meta']['nsessions']
        available_sessions = set(range(1, nsessions + 1))

        if required_sessions.issubset(available_sessions):
            # classify config from raw metadata
            if subset_dict[mouse]['conc']['b_meta']['intervals'][0][-1] == 1200:
                config_type.append(2)
            else:
                config_type.append(1)

            mslopes = []
            valid_mouse = True

            for session in ses_n:
                lick_df, _, _, bmeta, _ = sf.extract_data(data_dict, mouse, session)
                _, _, fit_parameters = matching_calc(lick_df, bmeta, session, plot_ind=False)

                if len(fit_parameters) > 0:
                    mslopes.append(fit_parameters[0][0])
                else:
                    valid_mouse = False
                    print(f'{mouse} skipped: no valid fit for session {session}')
                    break

            if valid_mouse and len(mslopes) == len(ses_n):
                all_slopes.append(mslopes)
            else:
                config_type.pop()

        else:
            print(f'{mouse} skipped: only has {nsessions} session(s)')

    ses1_sens = [slopes[0] for slopes in all_slopes]
    ses2_sens = [slopes[1] for slopes in all_slopes]

    if plot and len(all_slopes) > 0:
        plot_multises_slopes(ses1_sens, ses2_sens, config_type)

    return all_slopes, config_type

### PLOTTING 

def plot_multises_slopes(ses1_sens,ses2_sens,config_type):
    """
    Plot matching sensitivity across two sessions for multiple mice.

    Parameters
    ----------
    ses1_sens : array-like
        Sensitivity values for session 1.
    ses2_sens : array-like
        Sensitivity values for session 2.
    config_type : list of int
        Configuration label for each mouse, used to determine point outline
        color in the plot.

    Notes
    -----
    Each mouse is represented by a line connecting its sensitivity across the
    two sessions. Boxplots summarize the distribution in each session, and
    individual points are overlaid with outline colors indicating configuration
    type (e.g. config 1 vs config 2). This function is intended for visual
    comparison of within-mouse changes across sessions.
    """
    
    hfont = {'fontname':'Arial'} 
    fig,ax = plt.subplots(figsize=(3,3))

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    
    [sns.lineplot(x=[0,1],y=[ses1_sens[i],ses2_sens[i]],ax=ax,color='gainsboro',linewidth=2,zorder=-1) for i in range(0,len(ses1_sens))]
    
    dot_colors = ['silver' if c == 1 else 'orangered' for c in config_type]
    
    sns.boxplot(y=ses1_sens,x=0,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
    sns.boxplot(y=ses2_sens,x=0.2,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
    
    sns.stripplot(y=ses1_sens,x=0,edgecolor=dot_colors,facecolor='white',linewidth=2,zorder=-1)
    sns.stripplot(y=ses2_sens,x=0.2,edgecolor=dot_colors,facecolor='white',linewidth=2,zorder=-1)
    
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
