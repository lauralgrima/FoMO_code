import numpy as np
import pandas as pd

def model_data_df(multi_lick_df,video_df,session_n):
    '''Creates a dataframe with video data, reward schedule for each port, and visit time at each port at the sampling rate of the camera frame rate'''

    # get data just for one session
    lick_df1  = multi_lick_df.loc[multi_lick_df['session_number']==session_n]
    video_df1 = video_df.loc[video_df['session_number']==session_n]
    
#    framerate = 200 # might not always be the case, check based on number of frames
    
    
    ### EXTRACT USEFUL STUFF FROM LICK_DF AND PUT INTO VIDEO_DF - everything in video framerate timeframe 
    # dataframe of reward availability, port x event
    RA_updates_df = lick_df1[lick_df1['rewarded'].isnull()][['RA_1','RA_2','RA_3','RA_4','RA_5','RA_6','video_i']]
    RA_updates_df = RA_updates_df.drop(RA_updates_df[(RA_updates_df[['RA_1','RA_2','RA_3','RA_4','RA_5','RA_6']]==0).all(axis=1)].index)
    
    RA_updates_df.set_index('video_i',inplace=True)
    RA_updates_df = RA_updates_df.loc[RA_updates_df.index != 0]
    
    # dataframe of visits, port x visit  
    port_visits_df = lick_df1.loc[(lick_df1['unique_visit']==1)][['port','video_i']]    
    port_visits_long_df = pd.DataFrame(np.nan,columns=['port1','port2','port3','port4','port5','port6'],index=range(0,len(port_visits_df)))
    for i,port in enumerate(port_visits_df['port']):
        port_visits_long_df.iloc[i,int(port)-1] = 1
    port_visits_long_df['video_i'] = port_visits_df['video_i'].values
    port_visits_long_df.set_index('video_i',inplace=True)
    port_visits_long_df = port_visits_long_df.loc[port_visits_long_df.index != 0]
    
    # whether visit is rewarded (1 or 0)
    rewarded_df = lick_df1[['rewarded','video_i']]
    rewarded_df.set_index('video_i',inplace=True)
    rewarded_df = rewarded_df.loc[rewarded_df.index != 0]
    rewarded_df = rewarded_df.groupby(rewarded_df.index).first()
    
    # add to video_df for session
    video_df1 = pd.concat([video_df1,RA_updates_df,port_visits_long_df,rewarded_df],axis=1)

    return(video_df1)