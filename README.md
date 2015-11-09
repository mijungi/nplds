The core functions for running the VBEM algorithm over the class of models specified in the article are within the core_functions folder. We also use a toolbox for initializing the time scale of the GP kernel, that is in the folder gpml_matlab.

Within the Figure2_simulated_data, one can find functions to generate simulated data, then fit 3 different candidate models to it.

For figure 3 and 4 of the paper we fitted real neurophysiological data with the model. The data was gathered by Ecker et al. (State dependence of noise correlations in macaque primary visual cortex. Neuron, 82(1):235–248,2014)

The function figure3_4_main runs the fitting with multiple latent dimensionalities (k), and multiple different left-out trial sets for cross-validation. The procedure depends on multiple design and initialization parameters within files figure3_4_main, runAlexsdata_prediction_func and VBcomputeh_C_NSFR, these are marked via the "FITPARAM" comment. The directory paths to specify folders for input/output need to be changed to accomodate local environments, these are marked via "HARDCODED_DIRECTORY_PATHS" comments.

The outputs are going to be saved in the specified folders within a subfolder specifying the dimensionality and the run number. Firstly the initialization parameters are stored in init_data.mat, the updates of each iteration are in datastruct_MODELNAME_iter_n.mat, then the final fit are in final_MODELNAME_data.mat

Part of the analysis is evaluating the predictive performance on left-out trials. This is conducted by Figure4/analyze_predictions.m script, predicting parameters on left-out trials, then simulating a number of examples (set by nsamps within make_predictions.m) based on those params. We get the value of the error by comparing the mean (or mode) of the simulations to the real data. The error values will be stored for each dimensionality and run in matrix variables NSFR_error and PLDS_error.

For the illustration purposes, further analysis was performed mostly on the run with the best predictive performance in the dataset, which was used to create Figure 4/B,C.

The figures of the paper, barring visualization changes, can be reproduced exactly via running the create_figure3_and_4_final.m script file.

The code uses the GP-toolbox GPML version 3.6 by Carl Rasmussen and Chris Williams, http://www.gaussianprocess.org/gpml/code/matlab/doc/. 
