"""
Project: Diversification Rates and ADE
Description:
Calculate log odds of survival at a given age
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
from sklearn.linear_model import LogisticRegression

species = pd.read_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species_raw.txt", sep = "\t")

# assign an age to each species through the time bins it existed in, using 5 Myr time bins

species_lst = []
bin_max = []
bin_min = []
bin_mid = []
bin_no = []

for i, j in enumerate(species["species"]):
    if species["age"].values[i] < 5:
        species_lst.append(j)
        bin_max.append(int(math.ceil(species["fad"].values[i])))
        bin_min.append(int(math.ceil(species["lad"].values[i])))
        bin_mid.append(int(math.ceil(species["fad"].values[i])) - 2.5)
        bin_no.append(1)
    else:
        time_int_max = list(range(int(math.ceil(species["fad"].values[i])), int(math.ceil(species["lad"].values[i]) + 4), -5))
        time_int_min = list(range(int(math.ceil(species["fad"].values[i]) - 5), int(math.ceil(species["lad"].values[i] - 1)), -5))
        bin_num = 1
        for k in range(len(time_int_max)):
            species_lst.append(j)
            bin_max.append(time_int_max[k])
            bin_min.append(time_int_min[k])
            bin_mid.append(time_int_max[k] - 2.5)
            bin_no.append(bin_num)
            bin_num += 1

df_bins = pd.DataFrame()
df_bins["species"] = species_lst
df_bins["bin_max"] = bin_max
df_bins["bin_min"] = bin_min
df_bins["bin_mid"] = bin_mid
df_bins["age"] = bin_no

# assign survival (0) to all time bins apart from the last one in each species, which will be marked as extinction (1)

surv = []

for i, j in enumerate(df_bins["species"].values):
    if i < len(df_bins["species"].values)-1:
        if j == df_bins["species"].values[i + 1]:
            surv.append(0)
        else:
            surv.append(1)
    else:
        surv.append(1)

df_bins["survival"] = surv


time_int_max = list(range(145, 4, -5))

# calculate logistic regression between survival and age for each time bin

lgits = []
bin = []

l = 0 # the number of species doesn't add up - CHECK!
for i, j in enumerate(time_int_max):
    # if i < len(time_int_max) - 1:
    tb = df_bins.loc[df_bins["bin_max"] == j]
    l += len(tb)
    print(len(tb))
    if len(tb) > 1 and len(np.unique(tb["survival"].values)) > 1: # both 1 and 0 need to be present
        logreg = LogisticRegression()
        lgit = logreg.fit(np.array(tb["age"].values).reshape(-1, 1), np.array(tb["survival"].values).ravel())
        lgits.append(lgit.coef_[0][0])
        bin.append(j)

df_odds = pd.DataFrame()
df_odds["bin"] = bin
df_odds["log_reg"] = np.exp(lgits)
df_odds["log_log_reg"] = np.log(df_odds["log_reg"])

log_lgits = np.log(lgits)





