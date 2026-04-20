import h5_funcs as h5
from pathlib import Path

subject_IDs = ['6PG5','6PG6','6PG8']
regions     = ['NAc','DMS']
tasks       = ['conc']

def load_data(df_folder, subject_IDs, tasks):
    """
    Load behavioral, photometry, and video data for a set of subjects across
    one or more tasks.

    This function infers each subject's recording region from the filename of
    the corresponding `_behav` file rather than requiring regions to be passed
    manually. For each subject-task pair, it searches for a file matching:

        <subject_ID>_*_<task>_behav*

    Example filename:
        6PG33_DMS_conc_behav

    The region is extracted as the second underscore-delimited field in the
    filename (e.g., 'DMS').

    Parameters
    ----------
    df_folder : str or Path
        Path to the directory containing the data files.
    subject_IDs : list of str
        Subject IDs to load.
    tasks : list of str
        Task names to load, e.g. ['conc', 'prob'].

    Returns
    -------
    data_dict : dict
        Nested dictionary structured as:

            data_dict[subject_ID][task] = {
                'region'     : str,
                'mlick_df'   : DataFrame,
                'b_meta'     : dict or DataFrame,
                'mphoto_df'  : DataFrame,
                'photo_meta' : dict or DataFrame,
                'mvid_df'    : DataFrame,
                'v_meta'     : dict or DataFrame
            }

        Subject-task combinations with no matching `_behav` file are skipped.

    Notes
    -----
    - If multiple `_behav` files are found for a subject-task pair, the first
      match is used and a warning is printed.
    - If a filename does not match the expected convention, that entry is skipped.
    - This function assumes region is consistent within a subject-task dataset
      and can be inferred from the `_behav` filename.
    """
    df_path = Path(df_folder)
    data_dict = {}

    for subject_ID in subject_IDs:
        data_dict[subject_ID] = {}

        for task in tasks:
            behav_matches = list(df_path.glob(f"{subject_ID}_*_{task}_behav*"))

            if len(behav_matches) == 0:
                print(f"No _behav file found for {subject_ID}, task={task}")
                continue

            if len(behav_matches) > 1:
                print(f"Multiple _behav files found for {subject_ID}, task={task}:")
                for f in behav_matches:
                    print(f"  {f.name}")

            behav_file = behav_matches[0]
            parts = behav_file.stem.split("_")

            if len(parts) < 4:
                print(f"Unexpected filename format: {behav_file.name}")
                continue

            region = parts[1]

            multises_ldf, behav_metadata, photo_df, photo_metadata, video_df, video_metadata = h5.load_df(
                df_folder, subject_ID, region, task
            )

            data_dict[subject_ID][task] = {
                'region': region,
                'mlick_df': multises_ldf,
                'b_meta': behav_metadata,
                'mphoto_df': photo_df,
                'photo_meta': photo_metadata,
                'mvid_df': video_df,
                'v_meta': video_metadata
            }

    return data_dict