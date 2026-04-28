import numpy as np
import pandas as pd
import preprocessing.support_funcs as sf
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
from scipy import stats


def multi_travel(div_from_rand_filepath,data_dict,ses_n,plot=True):
    """
    Analyze travel behavior across mice for a selected session.
    
    Computes relationships between physical port distance, actual path length, and
    travel time, and compares normalized path length before versus after divergence
    from random behavior. Optionally plots correlation summaries and before/after
    path-length changes.
    
    Parameters
    ----------
    div_from_rand_filepath : str
        Path to divergence-from-random latency data.
    data_dict : dict
        Mouse/session behavioral and video data.
    ses_n : int
        Session number to analyze.
    plot : bool, default=True
        Whether to generate summary plots.
    """
        
    path_dists,pd_corrs = [],[]
    path_times,pt_corrs = [],[]
    dt_corrs = []
    mean_pathlen_rand,mean_pathlen_arand = [],[]

    diverge_lats = sf.import_div_time(filepath=div_from_rand_filepath)
    
    for i,mouse in enumerate(list(data_dict.keys())):
        print(mouse) # just to make sure running okay 
        
        # if just for one session vs. all sessions 
        if isinstance(ses_n,int): 
            tot_ses = data_dict[mouse]['b_meta']['nsessions'] # total number of sessions for given mouse 
            if tot_ses>=ses_n: # if session is in range for given mouse 
                lick_df,_,vid_df,bmeta,_ = sf.extract_data(data_dict,mouse,ses_n)  
                
                visits = lick_df[lick_df['unique_visit']==1]

                # distance from previous port to current port (physical, not travelled)  
                distances = [[0,14,18,70,72.2,65.5],[14,0,22.8,56,65.5,42],[18,22.8,0,72.2,70,56],
                             [70,56,72.2,0,18,22.8],[72.2,65.5,70,18,0,14],[65.5,42,56,22.8,14,0]]
                dist = pd.Series([distances[int(visits['port'].iloc[i])-1][int(visits['port'].iloc[i+1])-1] for i in range(0,len(visits)-1)]).shift(1)
                dist = dist.replace(0,1) # deal with zeros - can't divide by zero 
                
                # actual path length
                path_df = dist_trav_between_visits(lick_df,bmeta,vid_df)
                
                # normalise by distance
                norm_pathlens = path_df.iloc[1:-1]['path_length'].reset_index(drop=True)/dist[1:].reset_index(drop=True)

                # splitting for before and after  (for path length for the moment)
                diverge_lat = diverge_lats[diverge_lats['Mouse']==mouse]['Median latency (minutes)']
                visits      = visits.iloc[1:-1]
                if not np.isnan(diverge_lat.iloc[0]):
                    diverge_lat   = diverge_lat.iloc[0]*60 # in seconds 
                    
                    visits_rand  = visits[visits['event_time']<diverge_lat]
                    
                    path_rand    = norm_pathlens[:len(visits_rand)]
                    path_arand   = norm_pathlens[len(visits_rand):]

                    # calculate average path length
                    mean_pathlen_rand.append(path_rand.mean())
                    mean_pathlen_arand.append(path_arand.mean())

                path_dist_df,pd_corr = corr_dist_pathlen(path_df)
                path_time_df,pt_corr = corr_time_pathlen(lick_df,path_df)
                dist_time_df,dt_corr = corr_dist_time(lick_df,path_df)

        path_dists.append(path_dist_df)
        pd_corrs.append(pd_corr)
        
        path_times.append(path_time_df)
        pt_corrs.append(pt_corr)
        
        dt_corrs.append(dt_corr)
        
    if plot:
        plot_ind_dist(pd.concat(path_dists),True,'dist')
        plot_corr_coef([pd_corrs[i][0] for i in range(0,len(pd_corrs))])
        
        plot_ind_dist(pd.concat(path_times),True,'time') 
        plot_corr_coef([pt_corrs[i][0] for i in range(0,len(pt_corrs))])
        
        plot_corr_coef([dt_corrs[i][0] for i in range(0,len(dt_corrs))])
        
        plot_before_after(np.array(mean_pathlen_rand)/1000,np.array(mean_pathlen_arand)/1000,[]) # in metres



def multi_travel_multi_ses(data_dict,ses_n=[1,2],plot=True):
    """
    Compare mean travel path length across multiple sessions.
    
    For each mouse that has all requested sessions, computes the mean path
    length per session and optionally plots paired comparisons (e.g., session 1
    vs session 2).
    
    Parameters
    ----------
    data_dict : dict
        Mouse/session behavioral and video data.
    ses_n : list, default=[1, 2]
        Sessions to include in the comparison.
    plot : bool, default=True
        Whether to generate a before/after comparison plot.
    """
    all_pathlen = []
    for i, mouse in enumerate(list(data_dict.keys())):
        print(mouse)
        ses_count = range(1,data_dict[mouse]['b_meta']['nsessions']+1)
        ses_incl  = list(set(ses_count).intersection(ses_n)) # sessions to be included 
        if len(ses_incl)>1: # if animal has both sessions 1 + 2
            pathlen = []
            for session in ses_n:
                lick_df,_,vid_df,bmeta,_= sf.extract_data(data_dict,mouse,session)  
                
                # calculate
                path_df = dist_trav_between_visits(lick_df,bmeta,vid_df)
                pathlen.append(path_df['path_length'].mean())
            all_pathlen.append(pathlen)
            
    ses1_pathlen = [all_pathlen[i][0] for i in range(0,len(all_pathlen))]
    ses2_pathlen = [all_pathlen[i][1] for i in range(0,len(all_pathlen))]
            
    if plot: 
        plot_before_after(np.array(ses1_pathlen)/1000,np.array(ses2_pathlen)/1000)
                


### CALCULATING

def dist_trav_between_visits(lick_df,bmeta,vid_df):
    """
    Compute path length traveled between consecutive visits.
    
    Calculates the Euclidean distance traveled between successive unique visits
    using smoothed x/y coordinates from video data. The first visit has no
    preceding path and is assigned NaN.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing 'unique_visit', 'port', and 'video_i'.
    bmeta : dict-like
        Behavioral metadata containing at least 'animal_no'.
    vid_df : pandas.DataFrame
        Video data containing 'x_smooth' and 'y_smooth' coordinates.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
        - 'port': visited port
        - 'path_length': distance traveled to reach that port
        - 'animal': animal identifier
    """
    visits       = lick_df.loc[lick_df['unique_visit']==1][['port','video_i']]
    visits['video_i'].replace(0,np.nan,inplace=True) # change 0s to nans
    vid_visits   = vid_df.iloc[visits['video_i'].dropna()]
    start_frames = vid_visits.index
    end_frames   = np.array(vid_visits.index-1)[1:]
    start_frames = np.array(start_frames[:-1])
    
    distances    = []
    for i,frame in enumerate(start_frames):
        xy_df    = vid_df.iloc[start_frames[i]:end_frames[i]][['x_smooth','y_smooth']]
        distances.append(euclidean_dist(xy_df))

    path_df                = pd.DataFrame(visits['port'])
    distances              = np.insert(np.array(distances),0,np.nan)
    if len(path_df)>len(distances):
        distances = np.concatenate((distances,np.full((len(path_df)-len(distances-1)), np.nan)))
        
    path_df['path_length'] = distances
    path_df['animal']      = [bmeta['animal_no']]*len(path_df)
    
    return path_df 

def euclidean_dist(xy_df):
    """
    Compute total Euclidean distance of a trajectory.
    
    Calculates the cumulative distance traveled along a path defined by
    consecutive (x, y) coordinates.
    
    Parameters
    ----------
    xy_df : pandas.DataFrame or array-like
        Sequence of x and y coordinates (shape: n_samples × 2).
    
    Returns
    -------
    float
        Total Euclidean distance traveled along the trajectory.
    """
    total_distance = 0.0
    for i in range(1, len(xy_df)):
    # Compute Euclidean distance between consecutive points
        dist = np.linalg.norm(np.array(xy_df)[i] - np.array(xy_df)[i - 1])
        total_distance += dist

    return np.sum(np.sqrt(np.sum(np.diff(xy_df, axis=0)**2, axis=1)))
    
def transition_trajectories(lick_df,vid_df,port):
    """
    Extract transition trajectories into a specified port.
    
    For a target port, collects video trajectories from each preceding port into
    that target port, sorts trajectories by path length, and returns them grouped
    by preceding port.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing 'unique_visit', 'port', and 'video_i'.
    vid_df : pandas.DataFrame
        Video data containing 'x_smooth' and 'y_smooth' coordinates.
    port : int
        Target port for incoming transitions.
    
    Returns
    -------
    list
        List of length 6. Each element contains trajectory DataFrames for
        transitions from one preceding port into the target port, sorted from
        longest to shortest.
    """
    visits       = lick_df.loc[lick_df['unique_visit']==1][['port','video_i']].reset_index(drop=True)
    prior_vis    = visits.iloc[visits.loc[visits['port']==port].index-1]
    
    trajs_sorted = []
    for port in range(1,7): # will give trajectories to the function port from an individual port 
        prior_vidi   = prior_vis[prior_vis['port']==port]['video_i']
        curr_vidi    = visits.iloc[prior_vidi.index+1]['video_i']
        trajs        = [vid_df.iloc[int(prior_vidi.iloc[i]):int(curr_vidi.iloc[i])][['x_smooth','y_smooth']] for i, frame in enumerate(prior_vidi)]
        trajs = [traj for traj in trajs if not traj.empty]
        path_len     = [euclidean_dist(traj) for traj in trajs]
        trajs_sort   = sorted(zip(path_len,trajs),reverse=True) # longest to shortest 
        trajs_sorted.append([trajs_sort[i][1] for i in range(0,len(trajs_sort))])
        
    return trajs_sorted

def corr_dist_pathlen(path_df):
    """
    Compute correlation between port-to-port distance and path length.
    
    Adds a column of physical distances between consecutive ports and computes
    the Spearman correlation between this distance and the actual traveled path
    length.
    
    Parameters
    ----------
    path_df : pandas.DataFrame
        DataFrame containing at least:
        - 'port': visited port sequence
        - 'path_length': traveled distance between visits
    
    Returns
    -------
    tuple
        Updated DataFrame with 'dist' column, and Spearman correlation result
        (correlation coefficient and p-value).
    """

    ports     = path_df['port']
    distances = [[0,14,18,70,72.2,65.5],[14,0,22.8,56,65.5,42],[18,22.8,0,72.2,70,56],
                 [70,56,72.2,0,18,22.8],[72.2,65.5,70,18,0,14],[65.5,42,56,22.8,14,0]]
    
    dist = []
    for i in range(1,len(ports)):
        start_port = ports.iloc[i-1]
        end_port   = ports.iloc[i]
        dist.append(distances[int(start_port-1)][int(end_port-1)])
    
    dist            = np.insert(np.array(dist),0,np.nan)
    path_df['dist'] = dist
    
    corr = stats.spearmanr(path_df['dist'],path_df['path_length'],nan_policy='omit')

    return path_df,corr

def corr_time_pathlen(lick_df,path_df):
    """
    Compute correlation between time between visits and path length.
    
    Adds the time difference between consecutive unique visits and computes the
    Spearman correlation between travel time and path length.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing 'unique_visit' and 'event_time'.
    path_df : pandas.DataFrame
        DataFrame containing at least 'path_length' and corresponding visit rows.
    
    Returns
    -------
    tuple
        Updated DataFrame with 'time' column, and Spearman correlation result
        (correlation coefficient and p-value).
    """
    path_df['time'] = lick_df[lick_df['unique_visit']==1]['event_time'].diff()
    corr = stats.spearmanr(path_df['path_length'],path_df['time'],nan_policy='omit')
    
    return path_df,corr


def corr_dist_time(lick_df,path_df):
    """
    Compute correlation between port-to-port distance and time between visits.
    
    Adds physical distance between consecutive ports and time between visits,
    then computes the Spearman correlation between these variables.
    
    Parameters
    ----------
    lick_df : pandas.DataFrame
        Behavioral data containing 'unique_visit' and 'event_time'.
    path_df : pandas.DataFrame
        DataFrame containing at least 'port' and corresponding visit rows.
    
    Returns
    -------
    tuple
        Updated DataFrame with 'dist' and 'time' columns, and Spearman correlation
        result (correlation coefficient and p-value).
    """
    ports     = path_df['port']
    distances = [[0,14,18,70,72.2,65.5],[14,0,22.8,56,65.5,42],[18,22.8,0,72.2,70,56],
                 [70,56,72.2,0,18,22.8],[72.2,65.5,70,18,0,14],[65.5,42,56,22.8,14,0]]
    
    dist = []
    for i in range(1,len(ports)):
        start_port = ports.iloc[i-1]
        end_port   = ports.iloc[i]
        dist.append(distances[int(start_port-1)][int(end_port-1)])
    dist            = np.insert(np.array(dist),0,np.nan)
    path_df['dist'] = dist
        
    path_df['time'] = lick_df[lick_df['unique_visit']==1]['event_time'].diff()
    
    corr = stats.spearmanr(path_df['dist'],path_df['time'],nan_policy='omit')
    
    return path_df,corr


### PLOTTING

def plot_eg_traj(vid_df,time_lim=500,n_trajs=200,plot_colorbar=False):
    """
    Plot an example animal trajectory over a fixed time window.
    
    Plots smoothed x/y video coordinates, split into short trajectory segments
    with a time-varying colormap.
    
    Parameters
    ----------
    vid_df : pandas.DataFrame
        Video data containing 'x_smooth' and 'y_smooth' columns.
    time_lim : int, default=500
        Time window to plot, in seconds.
    n_trajs : int, default=200
        Number of samples per plotted trajectory segment.
    plot_colorbar : bool, default=False
        Whether to plot a horizontal colorbar for the trajectory colormap.
    """
    n_frames = time_lim*200
    
    plotx = vid_df['x_smooth'][:n_frames].dropna()
    ploty = vid_df['y_smooth'][:n_frames].dropna()
    
    if len(plotx)!=len(ploty):
        print('whoops')
        
    # divide into trajectories for getting plotting looking nice 
    sec_len = len(plotx)//n_trajs
    xsect   = [plotx[i *n_trajs:(i + 1) * n_trajs] for i in range(sec_len)]
    ysect   = [ploty[i *n_trajs:(i + 1) * n_trajs] for i in range(sec_len)]
    
    
    color   = cm.magma(np.linspace(0.9,0.1,num=1799))
    fig, ax = plt.subplots(figsize=(4.5,1.5))

    for trajectory in range(0,len(xsect)):
        ax.plot(xsect[trajectory],ysect[trajectory],color=color[trajectory],linewidth=2)
            
    ax.spines[['bottom','top','right','left']].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
            
    fig.tight_layout()
    
    if plot_colorbar:

        fig, ax = plt.subplots(figsize=(6, 1))

        norm = plt.Normalize(vmin=0.1, vmax=0.9)  # Adjust vmin and vmax as needed
        sm   = plt.cm.ScalarMappable(cmap='magma', norm=norm)
        sm.set_array([])  # You can provide an empty array-like object here
        fig.colorbar(sm, orientation='horizontal', ax=ax)
    


def plot_eg_transitions(trajs_sorted):
    """
    Plot example transition trajectories into a target port.
    
    For each preceding port, plots all trajectories leading into the target port,
    with one subplot per preceding port. Trajectories are colored for visual
    separation.
    
    Parameters
    ----------
    trajs_sorted : list
        Output from `transition_trajectories`. A list of length 6, where each
        element contains trajectory DataFrames (with 'x_smooth' and 'y_smooth')
        for transitions from a given port, sorted by length.
    """
    for prev_port in range(1,7):
        port_trajs = trajs_sorted[prev_port-1]
        
        color   = cm.magma(np.linspace(0.9,0.1,num=len(port_trajs)))
        fig, ax = plt.subplots(figsize=(4.5,1.5))

        for trajectory in range(0,len(port_trajs)):
            ax.plot(port_trajs[trajectory]['x_smooth'],port_trajs[trajectory]['y_smooth'],color=color[trajectory],alpha=0.5)
            
        ax.spines[['bottom','top','right','left']].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
            
        fig.tight_layout()

    
def plot_corr_coef(corrs_r):
    """
    Plot correlation coefficients across animals.
    
    Displays the distribution of correlation coefficients as a boxplot with
    overlaid individual data points.
    
    Parameters
    ----------
    corrs_r : array-like
        Correlation coefficients (e.g., one per animal).
    """
    fig,ax = plt.subplots(figsize=(3.5,3.5))
    hfont  = {'fontname':'Arial'} 
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    sns.boxplot(y=corrs_r,x=0,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
    sns.stripplot(y=corrs_r,x=0.1,edgecolor='black',facecolor='white',linewidth=2)
 
    ax.set_ylabel('corr. coeff. ($r$)',**hfont,fontsize=18)
    
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
    ax.set_yticklabels([0,0.2,0.4,0.6,0.8,1],fontsize=18)
    ax.set_xlim([-1,5])
    ax.set_ylim([0,1])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()


    
def plot_ind_dist(all_path_dists,no_return,measure):
    """
    Plot path length as a function of travel distance or travel time.
    
    Pools transitions across animals/sessions and visualizes either path length by
    physical port distance or path length by travel time.
    
    Parameters
    ----------
    all_path_dists : pandas.DataFrame
        Concatenated transition data containing 'port', 'path_length', and either
        'dist' or 'time'.
    no_return : bool
        Whether to exclude transitions where the animal returned to the same port.
    measure : {'dist', 'time'}
        Variable to plot against path length.
    """
    fig, ax = plt.subplots(figsize=(4.5,4))
    hfont  = {'fontname':'Arial'}
    
    # drop rows with nans
    all_path_dists = all_path_dists.dropna()
    
    # adding colour to all_path_dists 
    vcolor = cm.viridis(np.linspace(0.9,0.1,num=8))
 #   bl = [0.4,0.4,0.4,0.4]
    bl = [1,1,1,1]
    
    dot_colors = []
    for dist in all_path_dists['dist']:
        if dist == 0.0:
            dot_colors.append(bl)
        elif dist == float(14):
            dot_colors.append(vcolor[0])
        elif dist == float(18):
            dot_colors.append(vcolor[1])
        elif dist == float(22.8):
            dot_colors.append(vcolor[2])
        elif dist == float(42):
            dot_colors.append(vcolor[3])
        elif dist == float(56):
            dot_colors.append(vcolor[4])
        elif dist == float(65.5):
            dot_colors.append(vcolor[5])
        elif dist == float(70):
            dot_colors.append(vcolor[6])
        elif dist == float(72.2):
            dot_colors.append(vcolor[7])
            
    all_path_dists['dotc'] = dot_colors
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    if no_return:
        all_path_dists = all_path_dists[all_path_dists['port']!=all_path_dists['port'].shift(1)]

    if measure == 'time':
        
        ax.scatter(all_path_dists['time'],all_path_dists['path_length']/1000,alpha=0.5,color=all_path_dists['dotc'])
        ax.set_ylim([0.1,100])
        ax.set_yscale('log')
        ax.set_xlabel('travel time (s)',**hfont,fontsize=18)
        ax.set_ylabel('path length (m)',**hfont,fontsize=18)
        ax.set_yticks([0.1,1.0,10,100])
        ax.set_yticklabels([0.1,1.0,10,100],fontsize=18,**hfont)
        ax.set_xticks([1,10,100,1000])
        ax.set_xticklabels([1,10,100,1000],fontsize=18,**hfont)
        
        
    elif measure == 'dist':
        ax = sns.boxenplot(x=measure,y='path_length',data=all_path_dists[['path_length','dist']].dropna(),linewidth=1)
    
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(1.5)
        
    ax.yaxis.set_tick_params(width=1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')

    ax.spines[['left','bottom']].set_position(('outward', 10))

    fig.tight_layout()


def plot_before_after(before,after):
    """
    Plot paired comparisons between two conditions.
    
    Displays individual paired data points (connected by lines) along with
    boxplots and overlaid scatter points for two conditions (e.g., before vs
    after, or session 1 vs session 2).
    
    Parameters
    ----------
    before : array-like
        Values for the first condition (e.g., before).
    after : array-like
        Values for the second condition (e.g., after).
    """
    
    hfont = {'fontname':'Arial'} 
    fig,ax = plt.subplots(figsize=(3,3.5))

    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    
    [sns.lineplot(x=[1,3],y=[before[i],after[i]],ax=ax,color='gainsboro',linewidth=2) for i in range(0,len(before))]

    sns.boxplot(y=before,x=0,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
    sns.stripplot(y=before,x=0.1,edgecolor='black',facecolor='white',linewidth=2)
    
    sns.boxplot(y=after,x=0.2,ax=ax,color='black',fill=False,width=0.4,showfliers=False,showcaps=False)
    sns.stripplot(y=after,x=0.3,edgecolor='black',facecolor='white',linewidth=2)
    
    ax.set_ylabel('mean path length (m)',**hfont,fontsize=18)
    
    ax.set_yticks([0,1,2,3,4])
    ax.set_yticklabels([0,1,2,3,4],fontsize=18)
    ax.set_xlim([-1,5])
    ax.set_ylim([0,4])
    
    ax.spines[['bottom','top','right']].set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5)
    
    fig.tight_layout()

