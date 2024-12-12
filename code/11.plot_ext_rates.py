"""
Project: Diversification Rates and ADE
Author: Kristína Kocáková
Description:
Plot extinction rates calculated in Script 8
Plot Figure 4, S8
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

##########
#Figure 4
##########

#data calculated for species in a single bin (145 - 0.01)
rates_single = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_rates_single_bin.npy")
rates_single = rates_single[:, :, 5:850]
subset_bins_single = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_subset_bins_10_bins.npy")
subset_bins_single = subset_bins_single[5:850]
mean_rates_single = np.mean(rates_single, axis = 1)

rates_log_single = np.log(rates_single)

mean_cats_single = np.array([rates_log_single[0, :, 5], rates_log_single[0, :, 44], rates_log_single[0, :, 138], rates_log_single[0, :, 494]])
mean_cats_single = np.transpose(mean_cats_single)
mean_cats_single = np.mean(mean_cats_single, axis = 0)

#data calculated for 10 bins based on extinction regime

rates = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_rates_10_bins.npy")
ages = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_subset_bins_10_bins.npy")
mean_rates = np.mean(rates, axis = 1)
mean_rates = np.log(mean_rates)

#Extinction regimes
age_cats = [5, 44, 138, 494]

#from top to bottom
ext_events = [0, 2, 5, 6]
ext_episodes = [4, 8]
constant = [1, 3, 7, 9]

age_cat_ext = []
ext_reg = []
mean_rt_ext = []

for i in range(mean_rates.shape[0]):
    for j in age_cats:
        if i in ext_events:
            age_cat_ext.append(j)
            ext_reg.append("Extinction event")
            mean_rt_ext.append(mean_rates[i][j])
        elif i in ext_episodes:
            age_cat_ext.append(j)
            ext_reg.append("Extinction episode")
            mean_rt_ext.append(mean_rates[i][j])
        elif i in constant:
            age_cat_ext.append(j)
            ext_reg.append("Constant extinction period")
            mean_rt_ext.append(mean_rates[i][j])

dyct_ext = {"Age category": age_cat_ext, "Extinction regimes": ext_reg, "Mean extinction rate": mean_rt_ext}
df_ext = pd.DataFrame.from_dict(dyct_ext)


#ADE types

#ADE types
age_cats = [5, 44, 138, 494]

#from top to bottom
ADE0 = []
ADE1 = [9, 8, 4, 2, 1]
ADE2 = [7, 6, 5, 3, 0]

age_cat_ade = []
ade_cat = []
mean_rt_ade = []

for i in range(mean_rates.shape[0]):
    for j in age_cats:
    #     if i in ADE0:
    #         age_cat.append(j)
    #         ade_cat.append("ADE0")
    #         mean_rt.append(mean_rates[i][j])
        if i in ADE1:
            age_cat_ade.append(j)
            ade_cat.append("ADE1")
            mean_rt_ade.append(mean_rates[i][j])
        elif i in ADE2:
            age_cat_ade.append(j)
            ade_cat.append("ADE2")
            mean_rt_ade.append(mean_rates[i][j])

dyct_ade = {"Age category": age_cat_ade, "ADE category": ade_cat, "Mean extinction rate": mean_rt_ade}
df_ade = pd.DataFrame.from_dict(dyct_ade)


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
sns.set_style("white")

fig, axes = plt.subplots(3, 1, figsize = (6.69291/2, 8.26772))
sns.set_style("whitegrid")
# sns.set_context("notebook", font_scale=1)
fig.supxlabel("Longevity", fontsize=7)
fig.supylabel("Extinction Rate\n(log scaled)", fontsize=7)

for i in range(rates_single.shape[1]):
    axes[0].plot([0.5, 4.4, 13.8, 49.4], [rates_log_single[0, i, 5], rates_log_single[0, i, 44], rates_log_single[0, i, 138], rates_log_single[0, i, 494]], color = "#e3d5d5",  marker="o", markersize = 2.5, alpha = 0.07)
axes[0].plot([0.5, 4.4, 13.8, 49.4], mean_cats_single, color = "#260000",  marker="o", markersize = 2.5,) #"#260000"
ytcks = np.array(axes[0].get_yticks())
ytcks_vals = np.around(np.exp(ytcks), decimals = 2)
axes[0].set_yticks(ticks = list(ytcks), labels = list(ytcks_vals))
axes[0].set_xticks(ticks = [1], labels = [""])

sns.lineplot(ax = axes[1], data= df_ext, x = "Age category", y = "Mean extinction rate", hue = "Extinction regimes", palette={"Extinction event": "#260000", "Extinction episode": "#a81e1e", "Constant extinction period": "#f79797"}, marker= "o", linewidth = 3, markersize = 10, estimator='mean', errorbar='sd', err_style="bars")
ytcks = np.array(axes[1].get_yticks())
ytcks_vals = np.around(np.exp(ytcks), decimals = 2)
axes[1].set_yticks(ticks = list(ytcks), labels = list(ytcks_vals))
plt.xlabel("")
plt.ylabel("")
axes[1].set_xticks(ticks = [1], labels = [""])

sns.lineplot(ax = axes[2], data= df_ade, x = "Age category", y = "Mean extinction rate", hue = "ADE category", palette={"ADE1": "#FCBA03", "ADE2": "#DD614A"}, hue_order = ["ADE1", "ADE2"], marker= "o", linewidth = 3, markersize = 10, estimator='mean', errorbar='sd', err_style="bars")
ytcks = np.array(axes[2].get_yticks())
ytcks_vals = np.around(np.exp(ytcks), decimals = 2)
axes[2].set_yticks(ticks = list(ytcks), labels = list(ytcks_vals))
plt.xlabel("")
plt.ylabel("")
axes[2].set_xticks([5, 44, 138, 494], labels = ["Origin\n0.5 Ma", "Young\n4.4 Ma", "Middle aged\n13.8 Ma", "Elder\n49.4 Ma"])


sns.despine()
plt.tight_layout()


##########
#Figure S8
##########

rates = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_rates_single_bin.npy")
rates = rates[:, :, 5:850]
subset_bins = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_subset_bins_10_bins.npy")
subset_bins = subset_bins[5:850]
mean_rates = np.mean(rates, axis = 1)

rates_log = np.log(rates)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
sns.set_style("white")

fig, ax = plt.subplots(figsize = (6.69291, 8.26772/2))
sns.set_style("whitegrid")
# sns.set_context("notebook", font_scale=1)
fig.supxlabel("Longevity", fontsize=7)
fig.supylabel("Extinction Rate\n(log scaled)", fontsize=7)
for i in range(rates.shape[1]):
    ax.plot(subset_bins, rates_log[0, i], color = "#e3d5d5", alpha = 0.01) #"#260000"
ax.plot(subset_bins, np.mean(rates_log, axis=1)[0], color = "#260000")
ax.set_xticks(list(np.arange(0, 85, 10)))
ytcks = np.array(ax.get_yticks())
ytcks_vals = np.around(np.exp(ytcks), decimals = 2)
ax.set_yticks(ticks = list(ytcks), labels = list(ytcks_vals))