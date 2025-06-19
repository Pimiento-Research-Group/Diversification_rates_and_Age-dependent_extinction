"""
Project: Diversification Rates and ADE
Description:
Extraction of Weibull shape parameters predicted by ADE-Bayes, calculation of mean + 95% CIs
"""

from pandas import *
import numpy as np
import scipy.stats as st

data = read_csv("/path_to_pyrate_output_folder/ADE-Bayes/all_species_1_ADE_ADE_mcmc.log", sep = "\t")

np.mean(data["w_shape"])

st.norm.interval(0.95, loc=np.mean(data["w_shape"]), scale=st.sem(data["w_shape"]))

