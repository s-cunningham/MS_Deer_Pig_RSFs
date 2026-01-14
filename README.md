# MS Deer & Pig RSFs and Carrying Capacity

Code for resource selection functions for deer and pigs across the state of Mississippi. Data can be accessed via a public records request from Mississippi Department of Wildlife, Fisheries and Parks.

Manuscript currently under review and will be linked upon acceptance.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17575178.svg)](https://doi.org/10.5281/zenodo.17575178)

### Description of scripts:

`00_functions.R`

Helper functions

`01_setting_up_gps_data.R`

Reads in and combines multiple data files for deer and GPS data. Cleans locations, removing any that are outside of Mississippi. Exports the combined .csv file for each species.

`02_deer_home_ranges.R` and `02_pig_home_ranges.R`

There is one script for each species, but they do the same procedures.

Reads in combined data files and subsamples GPS tracts to 1 location every 4 hours. Contains hard-coded manual filtering of GPS data by date (to account for collars that fell off, etc.). Collars were split using Lavielle method in order to calculate an AKDE for each seasonal activity center. Some Lavielle splits were manually checked. Then, for each deer and segment, an AKDE was calculated. From the AKDEs, we built grids of 'available' locations, combining them into a single dataset along with 'used' locations. We scaled the weights that we will use in the RSF by individuals.

`03_deer_spatial_covariates.R` and `03_pigs_spatial_covariates.R`

For each used and available location, we extracted the proportion of different landcovers. The landcover values were saved alongside the 'used' and 'available' locations.

`04_deer_full_mixed_RSF.R` and `04_pigs_full_mixed_RSF.R`

These scripts contain the resource selection function models, implemented as fully mixed models (random intercepts & slopes) via `glmmTMB`. We extracted the regression coefficients for prediction.

`05_predicting_global_deer_map.R` and `05_predicting_global_pigs_map.R`

These scripts use the RSF coefficients to predict across the entire state of Mississippi. Due to the size of the landscape vs the resolution of the landscape rasters, we predicted across one county at a time and mosaicked them together. We also built rasters depicting the confidence intervals for the predictions. We sampled prediction rasters to 90 x 90 m.

`06_Fig1_landscape_map_with_gps.R`

Plotting NLCD values across Mississippi, with GPS locations on top.

`06_making_fig3.R`

Not actually Figure 3 anymore. Creates combined plot showing RSF maps, and boxplots of RSF coefficients.

`06_plotting_hab_patches.R`

Plotting patches as delineated according to levels of suitability by RAMAS.

`06_plotting_scaled_covariates.R`

`06_summarizing_regional_HSI_K.R`

Summarizing habitat suitability and carrying capacity based on RAMAS-generated patch types.

`07_read_RAMAS_output.R`

Reads in .txt files from RAMAS, summarizes per patch type. Plots of carrying capacity.

`08_deer_MSY_vital_rates.R` and `08_pigs_MSY_vital_rates.R`

Runs optimizer on literature demographic rates, constrained by estimated carrying capacity at two different densities.

`09_plotting_pop_projections.R`

Plots the population projections for the populations with optimized demographic rates.
