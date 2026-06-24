## Multisource detection, localization, and quantification (MDLQ) using continuous monitoring systems (CMS)

### Overview

This repository contains code to estimate methane emission source and rate on individual facilities using concentration data from point-in-space continuous monitoring systems. The code in this repository implements the MDLQ model described in "A Bayesian hierarchical model for methane emission source apportionment." The MDLQ model can accommodate scenarios in which multiple sources are emitting simultaneously. The code is separated into two main scripts located in the `code` directory: `code/MAIN_1_MDLQ.R` runs the MDLQ model and `code/MAIN_2_analyze_results.R` performs the analysis presented in the accompanying paper. The helper scripts contain auxiliary functions used to run the MDLQ model.

The following data are required to run the MDLQ model: 1) methane concentration measurements from a network of CMS sensors, and 2) simulated concentration measurements at the sensor locations from an atmospheric transport model. The latter must contain one set of simulated concentrations for each potential emission source on the site. The `input_data` directory contains the necessary data to run the MDLQ model on a subset of the ADED 2024 experiment conducted at the METEC facility in Fort Collins, Colorado. Specifically, the `input_data` directory contains four files. `input_data/sensor_locations.csv` and `input_data/source_locations.csv` contain the latitude, longitude, and height of the sensors and sources in the ADED 2024 experiment. This information is not directly needed to run the MDLQ model, but it is needed to forward simulate methane concentrations at the sensor locations. `input_data/input_data_SAMPLE.RData` contains methane concentration measurements from the CMS sensors and simulated methane concentrations at the sensor locations. Simulated concentrations were generated using the Gaussian puff atmospheric dispersion model as described in the accompanying manuscript and in Jia et al. (2025): https://doi.org/10.1038/s41598-025-99491-x. Code to run the Gaussian puff model can be found at: https://github.com/Hammerling-Research-Group/FastGaussianPuff. Note that due to the proprietary nature of the CMS sensor technology used in this study, only a sample of the CMS concentration data from the ADED 2024 experiment can be shared publicly; the `input_data/input_data_SAMPLE.RData` file contains only a week of CMS observations and the corresponding simulations. All of the CMS concentration measurements can be made available for research purposes upon request to the corresponding author. Finally, `input_data/ground_truth.RData` contains the true emission source and rate for the ADED 2024 experiment. This file is not necessary to run the MDLQ model, but it is necessary to perform the analysis in the `MAIN_2_analyze_results.R` file. It takes about 20 minutes to run the MDLQ model on the week-long sample data using 7 cores of an Apple MacBook Pro with an Apple M1 chip.

The `output_data` directory contains several pregenerated output files. Specifically, `output_data/MDLQ_output_30min_interval_30min_step_SAMPLE.RData` contains source and rate estimates from running the MDLQ model on the sample data provided in the `input_data` directory. `output_data/MDLQ_output_30min_interval_30min_step.RData` contains source and rate estimates from running the MDLQ model on the entire ADED 2024 dataset, which cannot be shared publicly as discussed above. Running the full output file through `code/MAIN_2_analyze_results.R` will reproduce the main results of the accompanying manuscript. The files starting with `output_data/pregenerated_` contain data objects used in `code/MAIN_2_analyze_results.R` that have been pregenerated to save time. If desired, these files can be reproduced by adjusting the appropriate toggles in `code/MAIN_2_analyze_results.R`. 


### Citation

If you use this code in your research, please cite the accompanying paper:

Daniels, Nychka, and Hammerling (2026). A Bayesian hierarchical model for methane emission source apportionment. *Annals of Applied Statistics*. DOI here.

BibTeX:
```
bibtex entry here.
```
