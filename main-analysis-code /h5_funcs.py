
# Taken from here: https://stackoverflow.com/questions/29129095/save-additional-attributes-in-pandas-dataframe/29130146#29130146

# Functions to store (h5store) and load (h5load) pandas dataframes in the hdf5 format. Particularly useful if working with large dataframes,
# or wanting to attribute metadata to dataframes.

import pandas as pd
import os

def h5store(filename, df, **kwargs):
    store = pd.HDFStore(filename)
    store.put('mydata', df)
    store.get_storer('mydata').attrs.metadata = kwargs
    store.close()

def h5load(store):
    data = store['mydata']
    metadata = store.get_storer('mydata').attrs.metadata
    return data, metadata


def load_df(df_folder,subject_ID,region,task):
    '''
    Loads dataframes and metadata saved in hdf5 format. Just for one mouse at a time. 
    
    Args: 
        root_folder          :   folder containing subfolders for behaviour, photometry, video data and dataframes
        subject_ID           :   string of subject whose data is to be loaded
        region               :   string, brain region to be loaded (e.g. 'NAc')
        task                 :   string, task to be loaded 

    Returns:
        A series of pairs of variables corresponding to dataframes and metadata. If the dataframe doesn't exist, will return empty list' 
    '''

    dfs = [dataframe for dataframe in os.listdir(df_folder) if not dataframe.startswith(".")] # ignore .DS_Store
    
    # assign variables as empty lists in case the dfs don't exist (behaviour should always be there)
    multises_df,behav_metadata,photo_df,photo_metadata,video_df,video_metadata = ([] for i in range(6))
    
    for df in dfs:
        
        df_subID,df_region,df_task = df[:-6].split("_")
        
        if df_subID == subject_ID and df_region == region and df_task == task:
             with pd.HDFStore(df_folder+'/'+subject_ID+'_'+region+'_'+task+df[-6:]) as store:
                 if df[-5:] == 'behav':
                     multises_df,behav_metadata = h5load(store)
                 elif df[-5:] == 'photo':
                     photo_df,photo_metadata = h5load(store)
                 elif df[-5:] == 'video':
                     video_df, video_metadata = h5load(store)

    return multises_df,behav_metadata,photo_df,photo_metadata,video_df,video_metadata 


