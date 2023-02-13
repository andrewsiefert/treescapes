## treescapes: Fitness landscapes for trees in the eastern US

This repository contains data and code to generate the results presented in: 

Siefert, A. and Laughlin, D.C. 2023. Estimating the net effect of functional traits on fitness across species and environments. *Methods in Ecology and Evolution* (in press).

1. Important files are organized in the `data`, `code`, and `results` folders (should be self explanatory).
2. Code to fit the demographic rate models (sapling survival, canopy tree survival, growth, and recruitment) are in the `code/demog` folder. The models are written in Stan (e.g., `growth_model.stan`) and fitted through R (e.g., `growth_model.R`) using the RStan package. The fitted models are exported as `.rds` files to the `results/models` folder. Be aware that the models take a long time to run (several days on a HPC cluster) and the fitted model objects are very large (many GB), so they aren't included in this repository, but can be reproduced using the data and code herein. 
3. The important parts of the fitted demographic rate models that are necessary for downstream analyses (i.e., posteriors of parameters of interest) are extracted from the fitted model objects using `code/demog/save_parameters_for_ipms.r` and saved in `results/demog_pars` folder (e.g., `growth_pars.rds`).
4. `code/ipm/performance_landscapes.r` uses the demographic rate models to calculate predicted demographic rates for various trait combinations at different temperatures and sizes, which are used to construct performance landscapes. 
5. `code/ipm/fitness_landscapes.r` uses the demographic rate models to construct IPMs and extract population growth rates (fitness) for variation trait combinations at different temperatures, which are used to construct fitness landscapes. 
6. `code/demographic_tradeoffs/demographic_tradeoffs.r` uses the demographic rate models to calculate demographic rates for trees with varying trait values to examine trait-mediated demographic trade-offs.
7. `code/model_checks` contains code to do posterior predictive checks and validation (e.g., Bayesian R2 and coverage) of demographic rate models using in-sample (training) and out-of-sample (test) data. 
