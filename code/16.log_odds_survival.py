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

# assign an age and survival (0 - survival, 1 - extinction) to each species through the time bins it existed in, using 5 Myr time bins

time_int_max = list(range(145, -1, -5))

species_lst = []
bin = []
age = []
surv = []

for j, k in enumerate(species["species"]):
    a = 1
    for i in time_int_max:
        if species["fad"].values[j] <= i and species["fad"].values[j] > i - 5:
            species_lst.append(k)
            bin.append(i)
            age.append(a)
            a += 1
            if species["lad"].values[j] < i - 5:
                surv.append(0)
            else:
                surv.append(1)
        elif species["fad"].values[j] >= i and species["lad"].values[j] <= i - 5:
            species_lst.append(k)
            bin.append(i)
            age.append(a)
            a += 1
            if species["lad"].values[j] < i - 5:
                surv.append(0)
            else:
                surv.append(1)
        elif species["lad"].values[j] <= i and species["lad"].values[j] >= i - 5:
            species_lst.append(k)
            bin.append(i)
            age.append(a)
            a += 1
            if species["lad"].values[j] < i - 5:
                surv.append(0)
            else:
                surv.append(1)

df_bins = pd.DataFrame()
df_bins["species"] = species_lst
df_bins["bin"] = bin
df_bins["age"] = age
df_bins["survival"] = surv

# calculate logistic regression between survival and age for each time bin

lgits = []
bin = []

for i in time_int_max:
    tb = df_bins.loc[df_bins["bin"] == i]
    print(len(tb))
    if len(tb) > 1 and len(np.unique(tb["survival"].values)) > 1: # both 1 and 0 need to be present
        logreg = LogisticRegression()
        lgit = logreg.fit(np.array(tb["age"].values).reshape(-1, 1), np.array(tb["survival"].values).ravel())
        lgits.append(lgit.coef_[0][0])
        bin.append(i)


df_odds = pd.DataFrame()
df_odds["bin"] = bin
df_odds["log_reg"] = lgits



sns.histplot(data = df_bins.loc[df_bins["bin"] == 55], x = "age", hue = "survival", multiple = "stack")
plt.show()

# plot these through time

