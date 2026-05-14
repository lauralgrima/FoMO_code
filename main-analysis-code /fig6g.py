import pandas as pd
import import_mat as im
import numpy as np
import matplotlib.pyplot as plt
from pypalettes import load_cmap


    
def MULTIalpha_KL(KL_path):
    '''
    Plot the relationship between learning rate (alpha) and KL divergence.

    Imports model-derived alpha and KL divergence values, concatenates data across
    mice, and visualizes the relationship between nominal learning rate (alpha)
    and divergence from the session 2 baseline. Lines connect values within mice,
    and points are colored by session.

    Parameters
    ----------
    KL_path : str
        Path to stored AQUA KL/alpha data.

    Returns
    -------
    None
        Generates a matplotlib figure showing alpha vs. KL divergence.
    '''
    _,alpha_KL    = im.import_KL_AQUA_data(KL_path)
    alpha_KL_conc = pd.concat(alpha_KL)
    alpha_KL_conc = alpha_KL_conc.iloc[:,:-5]
    
    ses_cols = [load_cmap("Starfish")(i) for i in [0.2,0.5,0.75]]
    
    hfont    = {'fontname':'Arial'} 
    fig, ax  = plt.subplots(figsize=(3.5,3)) 
    for i in range(0,6):
        alpha_byses = alpha_KL_conc[alpha_KL_conc.index==i][:3]
        KL_byses    = alpha_KL_conc[alpha_KL_conc.index==i][3:]
        [ax.plot(alpha_byses[j],KL_byses[j],color='gainsboro',alpha=0.5,linewidth=2,zorder=-1) for j in range(len(alpha_byses.iloc[0]))]
        [ax.scatter(alpha_byses[j],KL_byses[j],color=ses_cols,s=10) for j in range(len(alpha_byses.iloc[0]))]

    ax.set_xlim([-0.01,0.4])
    ax.set_xticks([0,0.4])    
    ax.set_ylim([0,0.07])
    ax.set_yticks([0,0.07])

    ax.set_xlabel(r'nominal $\alpha$',**hfont,fontsize='18')
    ax.set_ylabel('KL divergence from \n session 2 (s.d.)',**hfont,fontsize='18')
        
    ax.spines[['right','top']].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.spines[['left','bottom']].set_linewidth(1.5)
    
    fig.tight_layout()

def MULTIDA_KL(KL_path,incl_ses1=False):
    '''
    Plot the relationship between dopamine response and KL divergence across sessions.

    Imports paired dopamine (DA) response magnitudes and KL divergence values for
    defined session epochs, and visualizes their relationship for each mouse.
    Points are connected within mice and colored by session/epoch.

    Parameters
    ----------
    KL_path : str
        Path to stored AQUA DA/KL dataset. Expected format:
        columns correspond to DA responses followed by matched KL values
        (e.g., [end S1, start S2, end S2, start S3] and corresponding KL values).
    incl_ses1 : bool, optional
        If True, includes end-of-session-1 values in the analysis.
        If False, only uses session 2→3 transition epochs. Default is False.

    Returns
    -------
    None
        Generates a matplotlib figure showing DA response vs. KL divergence.
    '''
    DA_KL,_ = im.import_KL_AQUA_data(KL_path)
    
    hfont    = {'fontname':'Arial'} 
    fig, ax  = plt.subplots(figsize=(3.5,3)) 
    if incl_ses1:
        ses_cols = [load_cmap("Starfish")(i) for i in [0.25,0.5,0.75,1.0]]
    else:
       ses_cols = [load_cmap("Starfish")(i) for i in [0.2,0.5,0.75]]
       
    for i,mouse in enumerate(DA_KL):
        mDA_KL = DA_KL.iloc[mouse]
        if incl_ses1:
            xDA    = np.array(mDA_KL[:4])
            yKL    = np.array(mDA_KL[4:])
        else:
            xDA    = np.array(mDA_KL[1:4])
            yKL    = np.array(mDA_KL[5:])
        ax.plot(xDA,yKL,color='gainsboro',alpha=0.5,linewidth=2,zorder=-1)
        ax.scatter(xDA,yKL,color=ses_cols,s=20)
        
    ax.set_xlim([0,2.5])
    ax.set_xticks([0,2.5])    
    ax.set_ylim([-0.05,0.2])
    ax.set_yticks([0,0.2])
    
    ax.set_xlabel('session DA response',**hfont,fontsize='18')
    ax.set_ylabel('KL divergence from \n session 2 (s.d.)',**hfont,fontsize='18')
        
    ax.spines[['right','top']].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.tick_params(axis='both',which='both', labelsize=18,direction='in')
    ax.spines[['left','bottom']].set_linewidth(1.5)
    
    fig.tight_layout()

