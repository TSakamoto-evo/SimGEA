Code used in Sakamoto and Yeaman (2026). 
Each folder contains code for running different parts of the analysis.

The `stepping_stone` folder contains code for generating the simulation datasets used in this study.  
The `analyze_pvals` folder contains code for analyzing the GEA outputs.  
The `island` folder contains code for the island model analysis.  

`island/make_mig_matrix.py` generates the true migration rates.  
`island/coal_simu` contains the C++ simulator used to generate allele count data based on the true migration matrix.  
The output files were analyzed using both SimGEA and a method modified from Goel et al. (2026, bioRxiv; https://doi.org/10.64898/2026.01.19.700474).  
Code for running the Goel et al. (2026)-like analysis is stored in the `RoFi-like` folder.  
The Stan code was modified based on Goel et al. (2026)'s code available at 10.5281/zenodo.18307932.