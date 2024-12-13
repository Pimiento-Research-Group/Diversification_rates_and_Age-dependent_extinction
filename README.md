# 145-million-year global reconstruction of neoselachian diversification reveals new extinction events and persistent age-dependency 

This repository contains the code and data for the analyses presented in the manuscript:

**"145-million-year global reconstruction of neoselachian diversification reveals new extinction events and persistent age-dependency"**  
**Kristína Kocáková<sup>1*</sup>, Daniele Silvestro<sup>2,3,4*</sup>, Gregor Mathes<sup>1,5</sup>, Jaime A. Villafaña<sup>6</sup>, Catalina Pimiento<sup>*1,7</sup>**  

<sup>1</sup>Department of Paleontology, University of Zurich, Zurich, 8006, Switzerland<br />
<sup>2</sup>Department of Biology, University of Fribourg, Fribourg, 1700, Switzerland<br />
<sup>3</sup>Swiss Institute of Bioinformatics, Fribourg, 1700, Switzerland<br />
<sup>4</sup>Department of Biological and Environmental Sciences, Global Gothenburg Biodiversity Centre, University of Gothenburg, Gothenburg, 413 19, Sweden<br />
<sup>5</sup>GeoZentrum Nordbayern, Friedrich-Alexander University Erlangen-Nürnberg (FAU), Erlangen, 91054, Germany<br />
<sup>6</sup>Departamento de Ecología, Facultad de Ciencias, Universidad Católica de la Santísima Concepción, Concepción, 4090541, Chile<br />
<sup>7</sup>Department of Biosciences, Swansea University, Swansea, SA2 8PP, UK<br />

\*Corresponding authors: [kristina.kocakova@pim.uzh.ch](mailto:kristina.kocakova@pim.uzh.ch); [daniele.silvestro@unifr.ch](mailto:daniele.silvestro@unifr.ch); [catalina.pimientohernandez@pim.uzh.ch](mailto:catalina.pimientohernandez@pim.uzh.ch)

---

## Overview

This repository contains the scripts used to estimate rates of extinction and speciation, and evaluate age-dependent extinction (ADE) of neoselachian species and genera through the last 145 Myr.

### Repository Structure

- **`data/`**  
  Contains:
  - a .txt file with geological stage boundaries in milions of years (Myr)
  - PyRate input files generated and used in the analyses presented in the publication
  - four pre-trained models used for the estimation of ADE

- **`code/`**  
  Contains:
  - the scripts for the analyses
  - requirements.txt file with the required Python packages

---

## Usage notes

The PyRate program is required for a portion of the analyses in this study, the program and the instruction on how to compile it can be found in [this repository](https://github.com/dsilvestro/PyRate).
### Diversification rates
1. [**PyRate input step 1.**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/1.pyrate_input_step1.py)
  - filters the FINS dataset to produce a .txt file usable by the PyRate program, for an example of a PyRate input file see [this repository](https://github.com/dsilvestro/PyRate/wiki/1.-Preparing-input-file)
  - *input* - Occurrences data from the FINS Dataset openly accessible [here](https://zenodo.org/uploads/13983668)
  - *output* - .txt file formatted according to PyRate requirements
  - *alternative output* - if `adenn = True` argument is used, this script generates a .txt input file needed for the ADE-NN analyses
2. [**PyRate input step 2.**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/2.pyrate_input_step2.R)
  - converts the .txt file generated in Script 1. to a *_PyRate.py file readable by the PyRate program
  - the function in this script assigns an age to each occurrence randomy sampled from a uniform distribution delimited by the minimum and maximum age of the occurrence, for details see [this repository](https://github.com/dsilvestro/PyRate/wiki/1.-Preparing-input-file)
  - in our analyses we created an input file with 10 replicates of sampled ages
  - *input* - .txt file generated by Script 1.
  - *output* - *_PyRate.py file and a *_TaxonList.txt file containing a list of all unique taxa in the input file with the respective extinct/extant status. The taxon list file is not used in any further analyses here, but can be useful for data exploration. 
  - *Note*: As the process of the *_PyRate.py file creation has an element of stochasticity, we provide the files generated and used in subsequent analyses of this study in the [/data/ folder](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/tree/master/data) in this repository, for the purpose of reproducibility
3. **PyRate Analyses**<br />
- 3.1. [**Rate calculation**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/3.1.calculate_rates.sh)
    - estimates extinction, speciation and net diversification rates
    - the script is provided in a form of a shell script which can be submitted to a HPC Cluster, however it can also be ran locally by only running the final line of the script without the `srun` command and providing a replicate number for the `-j` argument
    - we used the following settings:
        - `-n` 50000000 - number of generations
        - `-s` 1000 - sample every 1000th generation
        - `-qShift` /path/to/ages.txt - allow peservation rate to vary in each geological stage (requires a .txt file, see *input* below)
        - `-min_dt` 0.5 - look for a rate shift every 0.5 Myr
        - `-mG`- allow preservation rate to vary amongst taxa
        - `-pP` 2 0 - prior for the gamma distribution of preservation rates
        - `-fast_burnin` 25000 -
        - `-fQ` 0.05 -
        - `-fU` 0.05 0.2 0 - 
        - `-singleton` 1 - excludes singletons, used in the main analysis and omitted in supplementary analysis
        - an overview of other PyRate settings can be accessed [here](https://github.com/dsilvestro/PyRate/wiki/2.-PyRate:-command-list#main-analysis-settings)
    - *input* - *_PyRate.py file generated in Script 2. and a .txt file containing starting and ending points of geological stages in Myr, which is provided [here](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/tree/master/data) as ages.txt, this file is required for the `-qShift` parameter
    - *output* - for each replicate produced in Script 2. the following files are generated - *_ex_rates.log, *_sp_rates.log, *_mcmc.log and *_sum.txt, here we had 10 replicates, so 40 output files in total are produced

- 3.2. [**Output processing**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/3.2.process_pyrate_output.sh)
  1. Combine the .log files into a single set of posterior values using the `-combLogRJ` command.
     - *input* - a path to the directory where the outputs from Script 3.1. are saved
     - *output* - 3 files - combined_10_mcmc.log, combined_10_ex_rates.log and combined_10_sp_rates.log
  2. Generate an .r script containing the mean rates and 95% CIs required for calculations and plotting using the `-plotRJ` command. This script also generates a default PyRate plot of the rates.
     - the following setting were used:
       - `-tag` combined - only use the combined posteriors, if not used a plot will be generated for each mcmc.log file in the directory
       - `-grid_plot` 0.1 - calculate the rate for each 0.1 Myr 
       - `-b` 0.1 - treat the first 10% of samples as burnin (ignore them)
     - *input* - a path to the directory where the outputs from Script 3.1. are saved or to the directory where the combined files are saved, if different
     - *output* - RTT_plots.r script and a .pdf plot
  3. Estimate speciation and extinction times using the `-ginput` command
     - *input* - combined_10_mcmc.log file
     - *output* - *_se_est.txt file containing the speciation and extinction times in Myr for each taxon in the dataset
- 3.3. [**Output cleaning**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/3.3.clean_pyrate_output.py)
  - this is an optional step created for a scenario in which the formatting of the output files isn't correct, which may happen when using an HPC Cluster (see notes directly in the script)
  - if required, this step should be done before step 3.2.
4. [**Plot the diversification rates**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/4.plot_rates.py)
  - plot diversification rates, produces figures 1., S1, S2 and S3
  - *input* - RTT_plots.r file generated in Script 3.2.
  - *output* - .pdf plots and an .xlsx file containing the rates at each 0.1 Myr
5. [**Find extinction and origination events**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/5.find_ext_events.py)
  - this script prints the time periods in Myrs where the extinction or origination rate reaches the intensity specified
  - the intesity of interest is defined in the `time_ext[rate_ext > 3*harm_mean]` line, here we looked for events of moderate intesity (i.e. at least 3 times higher than background), and intense intesity (i.e. at least 6 times higher than background, `time_ext[rate_ext > 6*harm_mean]`)
  - *input* - RTT_plots.r file generated in Script 3.2.
  - *output* - Myrs during which the desired level of rate intensity was reached, these will be printed in the console
---
### Age-dependent extinction (ADE)
6.  [**Train and evaluate ADE-NN models**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/6.model_training.py)
   - this script contains all the functions to simulate data, train, validate, test and save the ADE-NN models, and also to evaluate the performance of the models
   - **IMPORTANT**: The current version of the script produces the models in a Keras 2 format, which is not readable by Keras 3, therefore an environment with Keras 2 is needed. Alternative scripts producing models in a format compatible with Keras 3 will be added in the future.
   - *input* - none
   - *output* - 4 directories, each containing a trained model - a classifier and a predictive model for each of the three ADE types (ADE0, ADE1, ADE2 - see manuscript)
   - *Note*: This step can be skipped, pre-trained models are available [here](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/tree/master/data)
7. [**Predict ADE in empirical data**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/7.estimate_ADE.py)
  - estimate the class of ADE and the Weibull distribution shape parameters describing the distribution of longevities within a given time bin, including the prediction error. Sample size of each time bin is printed in the console during the estimation process
  - *inputs*:
    - path to the directories containing the trained models
        - pre-trained models used in this study are available [here](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/tree/master/data)
    - path to the ADE-NN input file - a .txt file generated using Script 1. and applying the `adenn = True` argument
    - a list defining the time bins of interest - this list can be based on extinction regimes defined in Script 5, or can represent any arbitrary time bins, such as geological epochs or stages
        - in this study we used the following time bins:
          - global species assemblage extinction regimes (results presented in Figure 3 and S4) - `time_slice = [[145, 95.610], [95.610, 83.502], [83.502, 73.096], [73.096, 71.995], [66.591, 65.791], [65.791, 55.785], [55.785, 38.774], [38.774, 33.471], [33.471, 3.452], [3.452, 0.01]]`
          - Cretaceous + Cenozoic (Figure 5a) - `time_slice = [[145, 0.01]]`
          - Geological periods (Figure 5b+c) - `time_slice = [[145, 66], [66, 23.03], [23.03, 0.01]]`
          - time bins used in Guinot & Condamine 2023 (Figure S6) - `time_slice = [[93.9, 66],[72.1, 66], [66, 56]]`
    - *output* - a 3D array containing the estimated values, saved as a .npy file, first dimension will reflect the number of time bins, second will reflect the number of replicates (we use 100), third is fixed (21 sets of values are estimated by the prediction function)
8. [**Estimate extinction rate as a function of age**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/8.calculate_ext_rate.py)
    - calculate extinction rates based on the Weibull shape parameter modelled in Script 7.
    - *input* - *.npy file from Script 7.
    - *output* - a 3D array of the estimated extinction rates per 0.1 Myr saved as a .npy file. The first two dimensions will be the same as the input file, the third will be the same as the number of subset bins
9. [**Find Myr representing age categories**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/9.find_age_categories.py)
    - age categories are used in the plotting scripts below
    - *input* - .txt file with estimated origination and extinction times from Script 3.2.3
    - *output* - ages in Myr representing young, middle-aged and elder taxa are printed in the console
10. [**Plot ADE results**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/10.plot_ADE_results.py)
    - the script contains 4 sections, each used to plot a figure (Figure 3, S4, S5, S6) - see notes inside the script
    - *inputs*:
        - .npy file containing results from Script 7.
        - .npy file containing results from Script 8.
        - list of time bins used when generating results of Script 7.
        - .xlsx file generated in Script 4.
    - *output* - figures
11. [**Plot extinction rates as a function of age**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/11.plot_ext_rates.py)
    - the script contains 2 sections, each used to plot a figure (Figure 4 and S8) - see notes inside the script
    - *input* - .npy file containing results from Script 8.
    - *output* - figures
12. [**Plot distribution of longevities**](https://github.com/Pimiento-Research-Group/Diversification_rates_and_Age-dependent_extinction/blob/master/code/12.plot_longevities.py)
    - plot Figure S7
    - *input* - *_se_est.txt file generated in Script 3.2.3
    - *output* - figures
13. d
14. d
15. d

















