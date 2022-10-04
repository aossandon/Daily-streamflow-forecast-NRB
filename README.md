# Daily-streamflow-forecast-NRB
This repository provides the dataset for the Narmada River Basin (NRB), R, and STAN scripts for implementing the Bayesian Hierarchical Model Combination (BHMC) framework proposed in _Ossandón et al. (2022)_. 
## Dataset for the NRB
This dataset contains the files with time series of potential covariates, daily peak monsoon (July-August) streamflow used to post-process daily VIC streamflow forecast across the Narmada River basin network (five gauges), India, for the period calibration (2003-2018). It also contains a file with basic information (longitude, latitude, and area) for the gauges considered here and observed data and covariates at the Handia gauge for the peak monsoon season 2021. The potential covariates comprise daily VIC forecasted (1- to 10-day lead time), simulated streamflow from each gauge, and 1, 2, 3, and 4-days accumulated spatial average observed precipitation from the area between the station gauges from 1 to 10-day lead times. The observed streamflow and gridded precipitation data were obtained from the India Water Resources Information System (India-WRIS) and the India Meteorology Department (IMD).
## Scripts
### R
- Library.R: Install and load all the packages required. Also, it contains all the personal functions created for this implementation. 
- Calibration_stan_2step_2covar_Diff_LeadTimes.R: Fit the best candidate BHMC for different lead times (1 to 10-day lead time). As output, it generates the file Par_bayesian_model_gamma_best_mod_kLT.rds with the posterior distribution of the model parameters. Where k indicates the lead time.
- Get_forecast_2021.R: Generate the ensemble streamflow forecast for the peak monsoon season 2021. Before running this script, the script Par_bayesian_model_gamma_best_mod_kLT.rds must be run for 1 to 5-days lead time. 

```{r, include=FALSE}
a=2
```
