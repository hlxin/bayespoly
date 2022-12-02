# bayespoly
Bayesian optimization code for Bayesian-Optimization-Assisted Discovery of Stereoselective Aluminum Complexes for Ring-Opening Polymerization of Racemic Lactide. 

The BO folder contains the scripts and data to run the Bayesian reaction optimization. BO/BO.py runs one iteration of Bayesian optimization using the training dataset BO/data/Al/experiment_index.csv (for all descriptors except OHE) or BO/data/Al/experiment_index_ohe.csv (for OHE as descriptor) together with descriptor values of the molecular fragments for each descriptor type in the same directory. BO/regression.py generates the results of 5-fold CV based regression performance test. BO/search.py generates the results of search efficiency test with random search as the benchmark. BO/edbo contains the re-coded scripts of edbo package (https://github.com/b-shields/edbo) for Bayesian optimization implementation.

The descriptors folder contains the scripts and molecular files to generate values of the 4 types of descriptors. descriptors/autoqchem contains the re-coded scripts of autoqchem package (https://github.com/doyle-lab-ucla/auto-qchem) for DFT descriptor extraction. descriptors/log_files contains all the Gaussian output files of molecular fragments with DFT descriptor information embedded.

The interpretation folder contains the scripts and data for model interpretation. interpretation/DFT_descriptors_86.csv contains the DFT descriptor values of our final Al-complex dataset. interpretation/whole_moleculars_Vbur.csv contains the calculated buried volume values of all whole molecules. interpretation/attibution_analysis.py generates the summary plot of SHAP analysis and the the multivariate linear regression plot based on fragmental DFT descriptors.    



