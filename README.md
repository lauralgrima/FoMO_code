This repository contains analysis code, task code, and model code for the manuscript:

**A global dopaminergic learning rate enables adaptive foraging across many options**

LL Grima, Y Guo, L Narayan, AM Hermundstad, JT Dudman

(DOI coming shortly)

The task-code folder contains scripts that specify the pyControl task definitions used in the interval and probabilistic schedule variants as outlined in the manuscript. For more information see: http://pycontrol.readthedocs.io

The main-analysis-code folder contains scripts with functions to generate manuscript figures. All analysis code is in Python 3, and all required packages should come by default with the Anaconda Python distribution. Refer to FoMO_main.py as the parent script which calls all other functions. 

To use the analysis code, download or clone this repository to obtain the main-analysis-code folder and the other_data folder. Download the file FoMO_dataset.zip from: https://doi.org/10.25378/janelia.31807939. The analysis code expects the folders with main-analysis-code, other_data, and the data to be in the same directory. 
