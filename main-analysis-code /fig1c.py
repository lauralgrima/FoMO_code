import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
                
def eg_event_traces(data_dict, mouse='6PG6', cutoff=30, axes_off=True):
    """
    Plot example visit events, reconstructed reward availability, and photometry
    traces for each port for a single mouse.

    Parameters
    ----------
    data_dict : dict
        Nested dictionary containing loaded session data.
    mouse : str, default='6PG6'
        Subject ID to plot.
    cutoff : float, default=30
        Maximum time to plot, in minutes.
    axes_off : bool, default=True
        If True, hide axes for all subplots.

    Returns
    -------
    None
        Creates one matplotlib figure per port.
    """

    lick_df  = data_dict[mouse]['conc']['mlick_df']
    photo_df = data_dict[mouse]['conc']['mphoto_df']
    visits   = lick_df.loc[lick_df['unique_visit'] == 1, ['port', 'event_time', 'rewarded']]

    for port in range(6):
        portn = f'RA_{port+1}'

        min_lick_df = lick_df.loc[
            (lick_df['unique_visit'] != 1) & (lick_df['rewarded'] != 1)
        ]

        RA      = min_lick_df.loc[:, [portn, 'event_time']].dropna()
        RA_time = RA.loc[RA[portn] == 1, 'event_time'] / 60
        rew_vis_time = visits.loc[
            (visits['port'] == port + 1) & (visits['rewarded'] == 1),
            'event_time'
        ] / 60
        un_vis_time = visits.loc[
            (visits['port'] == port + 1) & (visits['rewarded'] == 0),
            'event_time'
        ] / 60

        # limit to cutoff
        RA_time      = RA_time[RA_time < cutoff]
        rew_vis_time = rew_vis_time[rew_vis_time < cutoff]
        un_vis_time  = un_vis_time[un_vis_time < cutoff]

        fig, axs = plt.subplots(3, figsize=(10, 3))

        # visit raster
        axs[0].eventplot(
            [rew_vis_time, un_vis_time],
            linelengths=[0.5, 0.25],
            lineoffsets=-0.125,
            colors=['navy', 'firebrick'],
            linewidth=2
        )
        axs[0].set_xlim(0, cutoff)
        fig.tight_layout()

        # reward availability reconstruction
        bin_edges = np.arange(0, cutoff, 0.005)
        diff = pd.Series(
            np.histogram(RA_time, bins=bin_edges)[0] -
            np.histogram(rew_vis_time, bins=bin_edges)[0]
        )

        nonzero_diff = diff[diff != 0]
        repeat_idx   = nonzero_diff[nonzero_diff.eq(nonzero_diff.shift())].index
        diff.loc[repeat_idx] = 0

        avail = diff.replace(0, np.nan).ffill().fillna(0).replace(-1, 0)

        axs[1].plot(avail, color='black', linewidth=1.5)
        axs[1].set_xlim(0, len(avail))

        # photometry
        signal = photo_df['signal'].iloc[:(130 * 60 * cutoff)]
        axs[2].plot(signal, color='darkgrey', linewidth=1)
        axs[2].set_xlim(-2000, len(signal))

        if axes_off:
            for ax in axs.ravel():
                ax.set_axis_off()