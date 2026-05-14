import re
import pandas as pd 
import numpy as np
from scipy.stats import sem, t

def subset_mice(data_dict, task='conc', include_opto=False, opto_only=False, config=1, region=None):
    """
    Subset mice by task, opto status, spatial configuration, and region.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary of loaded data.
    task : str, default='conc'
        Task to include.
    include_opto : bool, default=False
        If False, exclude mice with IDs starting with '6PO'.
    opto_only : bool, default=False
        If True, include only opto mice (IDs starting with '6PO').
    config : {1, 2, None}, default=1
        Spatial configuration to include.
        - 1 → intervals[-1] != 1200
        - 2 → intervals[-1] == 1200
        - None → no config filtering
    region : {'NAc', 'DMS', None}, default=None
        Region to include. If None, no region filtering is applied.

    Returns
    -------
    subset_dict : dict
        Filtered dictionary with same structure as input.
    """

    def matches_config(mouse_data):
        intervals = mouse_data[task]['b_meta']['intervals'][0]
        is_config2 = intervals[-1] == 1200

        if config == 1:
            return not is_config2
        elif config == 2:
            return is_config2
        else:
            return True

    def matches_region(mouse_data):
        if region is None:
            return True
        return mouse_data[task]['b_meta']['region'] == region

    subset_dict = {
        subj: {task: data[task]}
        for subj, data in data_dict.items()
        if (
            task in data
            and (
                subj.startswith('6PO') if opto_only
                else (include_opto or not subj.startswith('6PO'))
            )
            and matches_config(data)
            and matches_region(data)
        )
    }

    return subset_dict


def extract_data(data_dict, mouse, n_session):
    """
    Assign session-specific data from `data_dict` to variables for analysis.
    Missing photometry or video data are returned as empty lists.
    """
    conc = data_dict[mouse]['conc']

    # licking / behavior
    mlick_df = conc['mlick_df']
    lick_df = mlick_df.loc[mlick_df['session_number'] == n_session]

    # photometry (optional)
    photo_df = conc.get('mphoto_df', [])
    if isinstance(photo_df, pd.DataFrame):
        if 'session_number' in photo_df.columns:
            ses_photo_df = photo_df.loc[photo_df['session_number'] == n_session]
        elif n_session == 1:
            ses_photo_df = photo_df
        else:
            ses_photo_df = []
    else:
        ses_photo_df = []

    photo_meta = conc.get('photo_meta', [])

    # behavior metadata
    behav_meta = conc['b_meta']

    # video (optional)
    avid_df = conc.get('mvid_df', [])
    if isinstance(avid_df, pd.DataFrame):
        if 'session_number' in avid_df.columns:
            vid_df = avid_df.loc[avid_df['session_number'] == n_session]
        elif n_session == 1:
            vid_df = avid_df
        else:
            vid_df = []
    else:
        vid_df = []

    return lick_df, ses_photo_df, vid_df, behav_meta, photo_meta


def longform(data_list):
    '''
    Takes data and redistributes to one sample/second. 
    '''
    longform_data,_ = np.histogram(data_list,bins=range(10800+1))
    return longform_data

def sort_data_by_interval(data,intervals):
    '''
    Data should be a list of lists, length 6. 
    Sorts order of data by increasing interval duration.
    '''
    sorted_data = sorted(zip(intervals,data),key=lambda x: x[0])
    data_only   = [sorted_data[i][1] for i, data_list in enumerate(sorted_data)]
    return(data_only)


def port_2_rank(bmeta, n_ses):
    '''
    Returns a dictionary that allows conversion from port number per se to port rank,
    sorted by rank (1 → 6).
    '''
    if bmeta['intervals'][n_ses-1][0] % 1 == 0:
        interval_dict = {30:1, 60:2, 240:3, 1200:5, 2400:6}
    elif bmeta['intervals'][n_ses-1][0] == 0.59:
        interval_dict = {0.59:1, 0.5:2, 0.37:3, 0.36:4, 0.16:6}
    else:
        print('intervals not included here - double check')
        return None

    ses_intervals = bmeta['intervals'][n_ses-1]
    ranks = np.array([interval_dict[interval] for interval in ses_intervals])

    if bmeta['intervals'][n_ses-1][0] % 1 == 0:
        ranks[np.where(ranks == 3)[0][0]] += 1
    else:
        ranks[np.where(ranks == 6)[0][0]] -= 1

    port_rank = dict(zip(range(1, 7), ranks))

    # sort by rank value
    return dict(sorted(port_rank.items(), key=lambda x: x[1]))


def _get_file_info(filename):
    '''Given a filename return the subject ID, brain region, task, and date.'''
    [subject_ID, region, task, date, time] = filename.split("_")
    time = time[:5]
    return subject_ID, region, task, date, time      
        

def _str_get(s,start='^',end='$'):
    'Return chunk of string s between specified start and end sequence.'
    return re.search(start+'(.*?)'+end, s).group(1)       


def NormalizeData(data):
    '''
    Normalise between zero and 1 (i.e. min-max normalisation.)
    Data should be a pandas Series (should make this more flexible in future).
    '''
    #return (data - np.min(data)) / (np.max(data) - np.min(data))
    return (data - data.min()) / (data.max() - data.min())


def NormalizeNeg(data):
    '''
    Normalise between -1 and 1.
    currently accepts np.arrays...
    '''
    
    return [2*((data[i]-np.min(data))/(np.max(data)-np.min(data)))-1 for i in range(0,len(data))]

def standardise(data):
    '''
    Mean centres at zero and with a variance of 1.
    '''
    mean_data    = np.mean(data)
    std_dev_data = np.std(data)
    return (data-mean_data)/std_dev_data
    

def sliding_window(win_size,data_list):
    '''
    Generates a list of arrays of sliding windows through a given list of data.
    '''
    windowed_data = [(data_list[i:i+win_size]) for i in range(len(data_list) - win_size+1)]
    return windowed_data


def moving_window(win_len,data_series):
    '''
    Generates a list of arrays of moving windows through a given pandas Series of data.
    '''
    data_series   = data_series.reset_index(drop=True)
    windowed_data = [data_series[(i*win_len):((i+1)*win_len)] for i in range(0,int(np.round(len(data_series)/win_len)-1))]
    final_window  = data_series.iloc[windowed_data[-1].index[-1]+1:data_series.index[-1]] # last window that isn't of same size as others 
    windowed_data.append(final_window)
    
    if len(windowed_data[-1]) == win_len+1:
        windowed_data[-1] = windowed_data[-1].iloc[:-1]
    
    return windowed_data
    

def convolve(win_size,data_list,mode='valid'):
    '''
    Similar to sliding_window, but then calculates the mean (convolves). 
    Set mode to:
        valid    :    convolution product only given for points where signals overlap completely
        full     :    at end points of convolution signals do not overlap completely. So boundary effects may be seen 
        same     :    will match the length of data_list
    '''
    return np.convolve(data_list,np.ones(win_size)/win_size,mode=mode)


def confidence_interval(data_list,confidence=0.95):
    '''
    data_list needs to be a list of arrays. Length of list is number of runs, and length of each array is number of windows
    '''
    
    means     = []
    intervals = []
     
    for i in range(len(data_list[0])):
        win_av = [run[i] for run in data_list] # averages across all runs for a single window 
        
        a = 1.0 *np.array(win_av)
        n = len(a)
        m, se = np.mean(a), sem(a)
        i = se*t.ppf((1+confidence)/2.,n-1) 
        
        means.append(m)
        intervals.append(i)
        
    return np.array(means),np.array(means)-np.array(intervals),np.array(means)+np.array(intervals)
    
def search_sequence_numpy(arr,seq):
    """ Find sequence in an array using NumPy only.

    Parameters
    ----------    
    arr    : input 1D array
    seq    : input 1D array

    Output
    ------    
    Output : 1D Array of indices in the input array that satisfy the 
    matching of input sequence in the input array.
    In case of no match, an empty list is returned.
    """

    # Store sizes of input array and sequence
    Na, Nseq = arr.size, seq.size

    # Range of sequence
    r_seq = np.arange(Nseq)

    # Create a 2D array of sliding indices across the entire length of input array.
    # Match up with the input sequence & get the matching starting indices.
    M = (arr[np.arange(Na-Nseq+1)[:,None] + r_seq] == seq).all(1)

    # Get the range of those indices as final output
    if M.any() >0:
        return np.where(np.convolve(M,np.ones((Nseq),dtype=int))>0)[0]
    else:
        return []         # No match found


# photometry-oriented

def extract_window(event_sample_numbers,signal,sample_rate,tw_start=-2,tw_length=5):
    # tw_start and tw_length are input as seconds, so this needs to be converted to samples - photometry data collected at rate of 130 samples/s
    tw_start   = tw_start*sample_rate
    tw_length  = tw_length*sample_rate
    
    if isinstance(event_sample_numbers, pd.DataFrame):
        event_sample_numbers = event_sample_numbers.astype(int).to_numpy().flatten()

    all_idx,ignored_idx = [],[]
    for i,s in enumerate(event_sample_numbers):
        if s + tw_start >= 0 and s + tw_start + tw_length < len(signal): # checking to see if the whole window is in range
            idx_start = s + tw_start
            all_idx.append([int(idx_start),int(idx_start+tw_length)])
        else: 
            ignored_idx.append(i)

    return all_idx,ignored_idx




