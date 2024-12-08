"""
Project: SHARK-XT Rates and ADE
Author: Kristína Kocáková
Description:
Script 4.
Calculating extinction rates in pre-determined time bins based on the ADE mode
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

def calc_scale(shape, longevity):
    scale = np.array(longevity/(gamma(1 + 1/shape)))
    return scale

def get_prob_ext(h, bins, b=0.1):
    d_inst = h[np.where((bins >= (b - bin_size)) & (bins <= b))[0][0]]
    conditional_d = d_inst / np.sum(h[bins[1:] > b])
    return conditional_d

def calc_ext_prob(i, data):
    shapes_0 = data[i, :, 18]
    shapes_1 = data[i, :, 3]
    shapes_2_1 = data[i, :, 9]
    shapes_2_2 =  data[i, :, 12]
    longevities_0 = data[i, :, 2]
    longevities_1 = data[i, :, 6]
    longevities_2 = data[i, :, 15]
    scales_0 = calc_scale(shapes_0, longevities_0)
    scales_1 = calc_scale(shapes_1, longevities_1)
    scales_2_1 = calc_scale(shapes_2_1, longevities_2)
    scales_2_2 = calc_scale(shapes_2_2, longevities_2)
    clas = data[i, 0, 1]


    if clas == 2:
        shapes_x = np.array([shapes_2_1, shapes_2_2])
        scales_x = np.array([scales_2_1, scales_2_2])
        d = np.random.weibull(shapes_x.flatten(), (n_samples, len(shapes_x.flatten()))) * scales_x.flatten()

        # remove extreme outliers (if any)
        d = d.flatten()
        d = d[d < 100]
    elif clas == 1:
        shapes_x = shapes_1
        d = np.random.weibull(shapes_x, (n_samples, len(shapes_x))) * scales_1
        # remove extreme outliers (if any)
        d = d.flatten()
        d = d[d < 100]
    elif clas == 0:
        shapes_x = shapes_0
        d = np.random.weibull(shapes_x, (n_samples, len(shapes_x))) * scales_0
        # remove extreme outliers (if any)
        d = d.flatten()
        d = d[d < 100]

    bin_size = 0.1
    bins = np.linspace(0, 100, int(100 / bin_size))
    h = np.histogram(d, bins=bins, density=True)[0]

    subset_bins = bins[bins < 93] #max longevity in our data is 93.something
    rates = np.array([get_prob_ext(h, bins, b) for b in subset_bins])
    rates = rates / bin_size
    # rates_list.append(rates)
    # rates_arr= np.array(rates_list)
    return rates, subset_bins

def calc_rates_list(i, data):
    rates_list = []
    for j in range(100):
        x, subset_bins = calc_ext_prob(i, data)
        rates_list.append(x)
    rates_arr = np.array(rates_list)
    # rates_log = np.log(rates_arr)
    return rates_arr, subset_bins


# Determine the age categories from empirical data
# Only needs to be done once when working with a new dataset

#calculate quantiles of ages globally

ages_f = pd.read_csv("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/all_species_15Myr/combined/combined_10_se_est_species_names.txt", sep = "\t")

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
age_50 = age_s[443:541] # middle 50%

age_1_med = np.median(age_1)
age_50_med = np.median(age_50)
age_99_med = np.median(age_99)

# End of age category determination


# Calculate extinction rates

n_samples = 10000
bin_size = 0.1

data = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/ADE_NN/empirical/bins_version2/2024_new/all_species_Jun_2024.npy")
# data[7, :, 1] = 0
rates = np.zeros((data.shape[0], 100, 930))
j = 0
for i in reversed(range(data.shape[0])):
     x, subset_bins = calc_rates_list(i, data)
     rates[j] = x
     j += 1

f = os.path.join("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/ADE_NN/empirical/bins_version2/2024_new/rates_Jun_2024.npy")
np.save(file=f, arr=rates)
f1 = os.path.join("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/ADE_NN/empirical/bins_version2/2024_new/subset_bins.npy")
np.save(file=f1, arr=subset_bins)

# End of extinction rate calculation

# Plotting

################
### Figure 4 ###
################

#for all data in 1 time bin
rates = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_rates_single_bin.npy")
rates = rates[:, :, 5:850]
subset_bins = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_subset_bins_10_bins.npy")
subset_bins = subset_bins[5:850]
mean_rates = np.mean(rates, axis = 1)

rates_log = np.log(rates)

mean_cats = np.array([rates_log[0, :, 5], rates_log[0, :, 44], rates_log[0, :, 138], rates_log[0, :, 494]])
mean_cats = np.transpose(mean_cats)
mean_cats = np.mean(mean_cats, axis = 0)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
sns.set_style("white")

fig, ax = plt.subplots(2, 1, figsize = (6.69291, 8.26772))
sns.set_style("ticks")
# sns.set_context("notebook", font_scale=1)
fig.supxlabel("Longevity", fontsize=7)
fig.supylabel("Extinction Rate\n(log scaled)", fontsize=7)
for i in range(rates.shape[1]):
    ax[0].plot(subset_bins, rates_log[0, i], color = "#e3d5d5", alpha = 0.01) #"#260000"
ax[0].plot(subset_bins, np.mean(rates_log, axis=1)[0], color = "#260000")
ax[0].set_xticks(list(np.arange(0, 85, 10)))
ytcks = np.array(ax[0].get_yticks())
ytcks_vals = np.around(np.exp(ytcks), decimals = 2)
ax[0].set_yticks(ticks = list(ytcks), labels = list(ytcks_vals))

for i in range(rates.shape[1]):
    ax[1].plot([0.5, 4.4, 13.8, 49.4], [rates_log[0, i, 5], rates_log[0, i, 44], rates_log[0, i, 138], rates_log[0, i, 494]], color = "#e3d5d5",  marker="o", markersize = 2.5, alpha = 0.07)
ax[1].plot([0.5, 4.4, 13.8, 49.4], mean_cats, color = "#260000",  marker="o", markersize = 2.5,) #"#260000"
ax[1].set_xticks([0.5, 4.4, 13.8, 49.4], labels = ["Origin\n0.5 Ma", "Young\n4.4 Ma", "Middle aged\n13.8 Ma", "Elder\n49.4 Ma"])
ytcks = np.array(ax[1].get_yticks())
ytcks_vals = np.around(np.exp(ytcks), decimals = 2)
ax[1].set_yticks(ticks = list(ytcks), labels = list(ytcks_vals))
sns.despine()
plt.tight_layout()
