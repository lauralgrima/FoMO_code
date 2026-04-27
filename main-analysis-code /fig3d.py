import import_mat as im
import support_funcs as sf
import numpy as np
import fig2a

def eg_cumu_visits_AQUA(data_dict,mat_GLM_path,mouse='6PG6',ses_n=1,plot=True):
    '''
    Calculate cumulative AQUA-predicted visits for each port for one mouse/session.

    Same idea as cumu_visits, but adapted for AQUA model output.

    Parameters:
        data_dict: dict
            Full behavioral data dictionary.

        mat_GLM_path: str
            Path to AQUA/GLM .mat data.

        mouse: str
            Mouse ID.

        ses_n: int
            Session number to analyze. Uses 1-indexing.

        plot: bool
            If True, plot cumulative visits.

    Returns:
        port_times: list
            Visit times for each port, sorted by interval.

        cumu_vis: list
            Cumulative visit counts for each port, sorted by interval.
    '''
    
    subset_dict_NAc = sf.subset_mice(data_dict,config=1, region='NAc')
    subset_dict_DMS = sf.subset_mice(data_dict,config=1, region='DMS')
    _,AQUA_beh_NAc  = im.import_GLMmat_data(mat_GLM_path,subset_dict_NAc,'NAc')
    _,AQUA_beh_DMS  = im.import_GLMmat_data(mat_GLM_path,subset_dict_DMS,'DMS')
    AQUA_beh        = {**AQUA_beh_DMS,**AQUA_beh_NAc}
    
    lick_df,photo_df,vid_df,bmeta,pmeta = sf.extract_data(data_dict,mouse,ses_n)
    
    mAQUA                = AQUA_beh[mouse]['AQUA_vis_invis']
    mAQUA['visit_times'] = lick_df[lick_df['unique_visit']==1]['event_time'].reset_index(drop=True)
    mAQUA                = mAQUA[mAQUA['session_number']==ses_n]
    
    total_vis  = mAQUA.iloc[:,:6].sum()
    cumu_vis   = [np.arange(1,total_vis[i]+1,1) for i in range(len(total_vis))]
    
    port_times = [mAQUA[mAQUA[i]==1]['visit_times'] for i in range(0,6)]

    port_times = sf.sort_data_by_interval(port_times,bmeta['intervals'][ses_n-1])
    cumu_vis   = sf.sort_data_by_interval(cumu_vis,bmeta['intervals'][ses_n-1])
    
    if plot:
        fig2a.plot_cumu_visits(port_times,cumu_vis)
