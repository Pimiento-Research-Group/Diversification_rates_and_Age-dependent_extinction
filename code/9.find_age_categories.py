"""
Project: Diversification Rates and ADE
Description:
Calculation of age categories from empirical data - determination which Myr value represents Young, Middle-aged and Elder taxa
Input file - .txt file generated using the -ginput function inscript 3.2.
"""

import glob, os, argparse, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import gamma
import pandas as pd
np.set_printoptions(suppress=True, precision=3)
import matplotlib.patheffects as pe
import itertools


# Determine the age categories from empirical data
# Only needs to be done once when working with a new dataset

#calculate quantiles of ages globally

ages_f = pd.read_csv("/path_to_pyrate_output_folder/combined_10_se_est_species_names.txt", sep = "\t")

ages_calc = []
for i in range(len(ages_f["species"])):
    x = ages_f["ts"][i] - ages_f["te"][i]
    ages_calc.append(x)
ages_f["age"] = ages_calc

ages_f = ages_f.loc[ages_f["te"] != 0]

ages = ages_f["age"]

age_s = sorted(ages) #arrange from youngest to oldest

age_1 = age_s[:98] # youngest 1%
age_99 = age_s[(984 - 98):]  # oldest 1%

age_1_med = np.median(age_1)
age_50_med = np.median(age_s)
age_99_med = np.median(age_99)

