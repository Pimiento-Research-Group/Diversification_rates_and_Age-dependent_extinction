"""
Project: SHARK-XT Rates and ADE
Author: Kristína Kocáková
Description:
Plotting of the distribution of predicted Weibull distribution parameters
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import matplotlib as mpl
np.set_printoptions(suppress=True, precision=3)  # prints floats, no scientific notation
import re
import matplotlib.patheffects as pe


###############
#Plot Figure 1.
###############
time_slice = [[145, 95.19], [95.19, 83.29], [83.29, 73.18], [73.18, 71.98], [66.48, 65.78], [65.78, 55.97], [55.97, 38.77], [38.77, 33.47], [33.47, 16.15], [16.15, 4.25], [4.25, 0.0117]]

# mean extinction rate plot prep
# only when using continuous rate
# rates = pd.read_excel("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/all_species_15Myr/rates.xlsx")
rates = pd.read_excel("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/fast_burnin_10_2024/species/rates.xlsx")


time_e = rates["Time_e"].abs()
rates["Time_e"] = time_e
time_e = np.array(rates["Time_e"].to_list())
rate_e = np.array(rates["Rate_e"].to_list())
mean_rates = []
for i in time_slice:
    x = rate_e[np.where((time_e <= i[0]) & (time_e >= i[1])),]
    mean_rates.append(np.mean(x))


position = []
width = []
value = [max(rate_e)+0.1]*len(time_slice)

for i in time_slice:
    w = i[0] - i[1]
    width.append(w)
    pos = i[0] - w/2
    position.append(pos)

# shape distribution prep

data_raw = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species.npy")

data = np.array(data_raw)
data = data[::-1, :, :]

classes = data[:, 0, 1].tolist()

ci_min = []
for i, j in enumerate(classes):
    if j == 2:
        x = []
        x1 = np.mean(data[i, :, 10])
        x.append([x1]*100)
        # x1 = data[i, :, 10]
        # x.append([np.min(x1)] * 100)
        x2 = np.mean(data[i, :, 13])
        x.append([x2]*100)
        # x1 = data[i, :, 13]
        # x.append([np.min(x1)] * 100)
        x = list(itertools.chain.from_iterable(x))
        ci_min.append(x)
    else:
        # ci_min.append([np.mean(data[i, :, 4])]*100)
        ci_min.append([np.min(data[i, :, 4])]*100)

ci_max = []
for i, j in enumerate(classes):
    if j == 2:
        x = []
        x1 = np.mean(data[i, :, 11])
        x.append([x1] * 100)
        # x1 = data[i, :, 11]
        # x.append([np.max(x1)] * 100)
        x2 = np.mean(data[i, :, 14])
        x.append([x2] * 100)
        # x2 = data[i, :, 14]
        # x.append([np.max(x2)] * 100)
        x = list(itertools.chain.from_iterable(x))
        ci_max.append(x)
    else:
        # ci_max.append([np.mean(data[i, :, 5])] * 100)
        ci_max.append([np.max(data[i, :, 5])] * 100)


# Import extinction rates

rates_subs = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_rates_10_bins.npy")
subset_bins = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-NN_10_2024/species_subset_bins_10_bins.npy")

mean_rates_x = np.mean(rates_subs, axis = 1)
mean_rates_x = np.log(mean_rates_x)

alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]
subsection = alphabet[1 : len(time_slice)+1]
subsection_2 = ["M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "Y"]

ext_events = [0, 2, 5, 6]
ext_episodes = [4, 8]
constant = [1, 3, 7, 9]
palette = []

for i in range(data.shape[0]):
    if i in ext_events:
        palette.append("#260000")
    elif i in ext_episodes:
        palette.append("#a81e1e")
    elif i in constant:
        palette.append("#f79797")

palette = list(reversed(palette))

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

fig, ax = plt.subplot_mosaic(
"""
ABM
ACN
ADO
AEP
AFQ
AGR
AHS
AIT
AJU
AKV
""",
    figsize=(6.69291, 8.26772), sharey= False,
    gridspec_kw={"width_ratios" : [1, 3, 2]},
    # per_subplot_kw = {"BCDE":{"sharex": "True"}} # - check after matplotlib update
    )

sns.set_style("whitegrid")
ax["A"].barh(position, height=width, width=value, color = palette, edgecolor = "none", alpha = 0.5)
ax["A"].plot(rate_e, time_e, color="#bf3e36", linewidth=2,
        path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
ax["A"].set_xlim([-0.1, max(rate_e)+0.1])
ax["A"].invert_xaxis()
ax["A"].invert_yaxis()
sns.despine()
ax["A"].set_ylabel("Time (Ma)", fontsize = 7)
ax["A"].set_xlabel("Extinction rate", fontsize = 7)

plt.tight_layout()

#add epochs

epoch_chrono = [122.5, 83, 61, 44.95, 28.465, 14.1815, 3.9565, 1.29585, 0.00585]
bar_vals = [-0.1] * 9
width_ep = [45, 34, 10, 22.1, 10.9, 17.7, 2.753, 2.5683, 0.0117]
epoch_chrono_names = ["Early Cretaceous", "Late Cretaceous", "Paleocene", "Eocene", "Oligocene", "Miocene", "Pli",
                      "Ple", ""]

ax["A"].barh(epoch_chrono, height=width_ep, width=bar_vals, edgecolor="black", color="white")
for i in range(len(epoch_chrono)):
    if i < 6:
        ax["A"].text(-0.05, epoch_chrono[i], epoch_chrono_names[i], fontsize=6.5, ha='center', va='center',
                     rotation="vertical")
    else:
        ax["A"].text(-0.05, epoch_chrono[i] + 0.4, epoch_chrono_names[i], fontsize=6.5, ha='center', va='center',
                     rotation="horizontal")

# plot shape distributions
for i in range(data.shape[0]):
    if i in ext_events:
        color = "#260000"
    elif i in ext_episodes:
        color = "#a81e1e"
    else:
        color = "#f79797"
    if classes[i] == 1:
        sns.histplot(ax = ax[subsection[i]], x = data[i, :, 3], color = "None", edgecolor = "None", kde=True)
        ax[subsection[i]].lines[0].set_color(color)
        ax[subsection[i]].plot([ci_min[i][0], ci_max[i][0]], [6, 6], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].axhline(y = 0, color = color, lw = 3)
        ax[subsection[i]].set_ylabel("")
        ax[subsection[i]].set_xlim(0.2, 3)
        ax[subsection[i]].set_xticks([0.1, 0.5, 1, 1.5, 2, 2.5, 3], labels = ["", "", "", "", "", "", ""])
        rates_x = [mean_rates_x[i, 5], mean_rates_x[i, 44], mean_rates_x[i, 138], mean_rates_x[i, 494]]
        # ax[subsection[i]].sharex(ax["I"])
        ax[subsection_2[i]].plot([5, 44, 138, 494], rates_x, marker="o", markersize = 2.5, c=color, linewidth=1,
                   path_effects=[pe.Stroke(linewidth=3, foreground='w', alpha=.7), pe.Normal()])
        ax[subsection_2[i]].spines.left.set_visible(False)
        ax[subsection_2[i]].spines.top.set_visible(False)
        ax[subsection_2[i]].yaxis.tick_right()
        ax[subsection_2[i]].set_ylabel("")
        ax[subsection_2[i]].set_xticks([5, 44, 138, 494], labels = ["","","",""])
        # ax[subsection_2[i]].sharey(ax["Y"])
        ax[subsection_2[i]].set_ylim([-4 ,-0.4])
        ytcks = np.array(ax[subsection_2[i]].get_yticks())
        ytcks_vals = np.around(np.exp(ytcks), decimals=2)
        ax[subsection_2[i]].set_yticks(ticks=list(ytcks), labels=list(ytcks_vals))
    elif classes[i] == 2:
        sns.histplot(ax=ax[subsection[i]], x=data[i, :, 9], color = "None", edgecolor = "None", kde=True)
        ax[subsection[i]].lines[0].set_color(color)
        sns.histplot(ax=ax[subsection[i]], x=data[i, :, 12], color = "None", edgecolor = "None", kde=True)
        ax[subsection[i]].lines[1].set_color(color)
        ax[subsection[i]].plot([ci_min[i][0], ci_max[i][0]], [6, 6], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].plot([ci_min[i][101], ci_max[i][101]], [9, 9], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].axhline(y=0, color = color, lw=3)
        ax[subsection[i]].set_ylabel("")
        # ax[subsection[i]].sharex(ax["I"])
        ax[subsection[i]].set_xlim(0.2, 3)
        ax[subsection[i]].set_xticks([0.1, 0.5, 1, 1.5, 2, 2.5, 3], labels = ["", "", "", "", "", "", ""])
        rates_x = [mean_rates_x[i, 5], mean_rates_x[i, 44], mean_rates_x[i, 138], mean_rates_x[i, 494]]
        ax[subsection_2[i]].plot([5, 44, 138, 494], rates_x, marker="o", markersize = 2.5, c=color, linewidth=1,
                   path_effects=[pe.Stroke(linewidth=3, foreground='w', alpha=.7), pe.Normal()])
        ax[subsection_2[i]].spines.left.set_visible(False)
        ax[subsection_2[i]].spines.top.set_visible(False)
        ax[subsection_2[i]].yaxis.tick_right()
        ax[subsection_2[i]].set_ylabel("")
        ax[subsection_2[i]].set_xticks([5, 44, 138, 494], labels = ["","","",""])
        # ax[subsection_2[i]].sharey(ax["Y"])
        ax[subsection_2[i]].set_ylim([-4, -0.4])
        ytcks = np.array(ax[subsection_2[i]].get_yticks())
        ytcks_vals = np.around(np.exp(ytcks), decimals=2)
        ax[subsection_2[i]].set_yticks(ticks=list(ytcks), labels=list(ytcks_vals))
    elif classes[i] == 0:
        sns.histplot(ax = ax[subsection[i]], x = data[i, :, 18], color = "white",edgecolor = "None", kde=True)
        ax[subsection[i]].lines[0].set_color(color)
        ax[subsection[i]].plot([ci_min[i][0], ci_max[i][0]], [6, 6], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].axhline(y = 0, color = color, lw = 3)
        ax[subsection[i]].set_ylabel("")
        # ax[subsection[i]].sharex(ax["B"])
        ax[subsection[i]].set_xlim(0.2, 3)
        ax[subsection[i]].set_xticks([0.1, 0.5, 1, 1.5, 2, 2.5, 3], labels = ["", "", "", "", "", "", ""])
        rates_x = [mean_rates_x[i, 5], mean_rates_x[i, 44], mean_rates_x[i, 138], mean_rates_x[i, 494]]
        # ax[subsection[i]].sharex(ax["I"])
        ax[subsection_2[i]].plot([5, 44, 138, 494], rates_x, marker="o", markersize = 2.5, c=color, linewidth=1,
                   path_effects=[pe.Stroke(linewidth=3, foreground='w', alpha=.7), pe.Normal()])
        ax[subsection_2[i]].spines.left.set_visible(False)
        ax[subsection_2[i]].spines.top.set_visible(False)
        ax[subsection_2[i]].yaxis.tick_right()
        ax[subsection_2[i]].set_ylabel("")
        ax[subsection_2[i]].set_xticks([5, 44, 138, 494], labels = ["","","",""])
        # ax[subsection_2[i]].sharey(ax["Y"])
        ax[subsection_2[i]].set_ylim([-4, -0.4])
        ytcks = np.array(ax[subsection_2[i]].get_yticks())
        ytcks_vals = np.around(np.exp(ytcks), decimals=2)
        ax[subsection_2[i]].set_yticks(ticks=list(ytcks), labels=list(ytcks_vals))

    if subsection[i] == subsection[-1]:
        ax[subsection[i]].set_xlabel("Estimated direction of age-dependency", fontsize=7)
        ax[subsection[i]].set_xticks([0.5, 1, 1.5, 2, 2.5, 3], labels = [0.5, 1, 1.5, 2, 2.5, 3])
    if subsection_2[i] == subsection_2[-1]:
        ax[subsection_2[i]].set_xlabel("Age category", fontsize=7)
        ax[subsection_2[i]].set_xticks([5, 50, 138, 494], labels = ["Origin", "Young", "Middle aged", "Elder"], fontsize = 6, ha = "left")
    # axes[i, j].set_titles("")
    ax[subsection[i]].set(yticks=[])
    # axes[i, j].set_ylabels()
    # axes[i, j].despine(left=True, bottom=True)
    ax[subsection[i]].axvline(x=1, linestyle="dashed", color="black")

sns.despine(left=True, bottom=True)
plt.tight_layout()

###############
#Plot Figure S4
###############
# extinction rate with highlighted bins + shape distributions

time_slice = [[145, 94.885], [94.885, 86.282], [86.282, 73.777], [73.777, 71.576], [67.075, 65.374], [65.374, 39.164], [39.164, 33.362], [33.362, 3.351], [3.351, 0.001]]
#extinction rare plot prep
rates = pd.read_excel("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/subsets/selachii.xlsx")

time_e = rates["Time_e"].abs()
rates["Time_e"] = time_e
time_e = np.array(rates["Time_e"].to_list())
rate_e = np.array(rates["Rate_e"].to_list())

position = []
width = []
value = [max(rate_e)+0.1]*len(time_slice)

for i in time_slice:
    w = i[0] - i[1]
    width.append(w)
    pos = i[0] - w/2
    position.append(pos)

# shape distribution prep

data_raw = np.load("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/ADE_NN/empirical/bins_version2/2024_new/subsets/selachii.npy")

data = np.array(data_raw)
data = data[::-1, :, :]

classes = data[:, 0, 1].tolist()

ci_min = []
for i, j in enumerate(classes):
    if j == 2:
        x = []
        x1 = np.mean(data[i, :, 10])
        x.append([x1]*100)
        # x1 = data[i, :, 10]
        # x.append([np.min(x1)] * 100)
        x2 = np.mean(data[i, :, 13])
        x.append([x2]*100)
        # x1 = data[i, :, 13]
        # x.append([np.min(x1)] * 100)
        x = list(itertools.chain.from_iterable(x))
        ci_min.append(x)
    else:
        # ci_min.append([np.mean(data[i, :, 4])]*100)
        ci_min.append([np.min(data[i, :, 4])]*100)

ci_max = []
for i, j in enumerate(classes):
    if j == 2:
        x = []
        x1 = np.mean(data[i, :, 11])
        x.append([x1] * 100)
        # x1 = data[i, :, 11]
        # x.append([np.max(x1)] * 100)
        x2 = np.mean(data[i, :, 14])
        x.append([x2] * 100)
        # x2 = data[i, :, 14]
        # x.append([np.max(x2)] * 100)
        x = list(itertools.chain.from_iterable(x))
        ci_max.append(x)
    else:
        # ci_max.append([np.mean(data[i, :, 5])] * 100)
        ci_max.append([np.max(data[i, :, 5])] * 100)

alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
subsection = alphabet[1 : len(time_slice)+1]

# change depending on the subset
ext_events = [0, 2, 5, 6]
ext_episodes = [4, 8]
constant = [1, 3, 7, 9]
palette = []

for i in range(data.shape[0]):
    if i in ext_events:
        palette.append("#260000")
    elif i in ext_episodes:
        palette.append("#a81e1e")
    elif i in constant:
        palette.append("#f79797")

palette = list(reversed(palette))
"""
from palettable.cmocean.sequential import Amp_20
from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list("LaJolla_custom", Amp_20.mpl_colors[5:20])
palette = []
for i in mean_rates:
    palette.append(cmap(i*5))
"""

# start of plotting
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

fig, ax = plt.subplot_mosaic(
"""
AB
AC
AD
AE
AF
AG
AH
AI
AJ
AK
""",
    figsize=(6.69291, 8.26772), sharey= False,
    gridspec_kw={"width_ratios" : [1, 3]},
    # per_subplot_kw = {"BCDE":{"sharex": "True"}} # - check after matplotlib update
    )

sns.set_style("whitegrid")
ax["A"].barh(position, height=width, width=value, color = palette, edgecolor = "none", alpha = 0.5)
ax["A"].plot(rate_e, time_e, color="#bf3e36", linewidth=2,
        path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
ax["A"].set_xlim([-0.1, max(rate_e)+0.1])
ax["A"].set_yticks(ticks=[145, 100.5, 66, 56, 33.9, 23.03, 5.333, 0],
                   labels=[-145, -100.5, -66, -56, -33.9, -23.03, -5.333, 0])
ax["A"].invert_xaxis()
ax["A"].invert_yaxis()
sns.despine()
ax["A"].set_ylabel("Time (Ma)", fontsize = 9)
ax["A"].set_xlabel("Extinction rate", fontsize = 9)

plt.tight_layout()

#add epochs

epoch_chrono = [122.5, 83, 61, 44.95, 28.465, 14.1815, 3.9565, 1.29585, 0.00585]
bar_vals = [-0.1] * 9
width_ep = [45, 34, 10, 22.1, 10.9, 17.7, 2.753, 2.5683, 0.0117]
epoch_chrono_names = ["Early Cretaceous", "Late Cretaceous", "Palaeocene", "Eocene", "Oligocene", "Miocene", "Pli",
                      "Ple", ""]

ax["A"].barh(epoch_chrono, height=width_ep, width=bar_vals, edgecolor="black", color="white")
for i in range(len(epoch_chrono)):
    if i < 6:
        ax["A"].text(-0.05, epoch_chrono[i], epoch_chrono_names[i], fontsize=6.5, ha='center', va='center',
                     rotation="vertical")
    else:
        ax["A"].text(-0.05, epoch_chrono[i] + 0.4, epoch_chrono_names[i], fontsize=6.5, ha='center', va='center',
                     rotation="horizontal")

# plot shape distributions
#
for i in range(data.shape[0]):
    if i in ext_events:
        color = "#260000"
    elif i in ext_episodes:
        color = "#a81e1e"
    else:
        color = "#f79797"
    if classes[i] == 1:
        sns.histplot(ax = ax[subsection[i]], x = data[i, :, 3], color = "None", edgecolor = "None", kde=True)
        ax[subsection[i]].lines[0].set_color(color)
        ax[subsection[i]].plot([ci_min[i][0], ci_max[i][0]], [6, 6], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].axhline(y = 0, color = color, lw = 5)
        ax[subsection[i]].set_ylabel("")
        ax[subsection[i]].set_xticks(ticks=[],
                           labels=[])
        ax[subsection[i]].sharex(ax["B"])
    elif classes[i] == 2:
        sns.histplot(ax=ax[subsection[i]], x=data[i, :, 9], color = "None", edgecolor = "None", kde=True)
        ax[subsection[i]].lines[0].set_color(color)
        sns.histplot(ax=ax[subsection[i]], x=data[i, :, 12], color = "None", edgecolor = "None", kde=True)
        ax[subsection[i]].lines[1].set_color(color)
        ax[subsection[i]].plot([ci_min[i][0], ci_max[i][0]], [6, 6], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].plot([ci_min[i][101], ci_max[i][101]], [9, 9], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].axhline(y=0, color = color, lw=5)
        ax[subsection[i]].set_ylabel("")
        ax[subsection[i]].sharex(ax["B"])
        ax[subsection[i]].set_xticks(ticks=[],
                           labels=[])
    elif classes[i] == 0:
        sns.histplot(ax = ax[subsection[i]], x = data[i, :, 18],color = "None", edgecolor = "None",kde=True)
        ax[subsection[i]].lines[0].set_color(color)
        ax[subsection[i]].plot([ci_min[i][0], ci_max[i][0]], [6, 6], "--|", color="#888a89", alpha = .5)
        ax[subsection[i]].axhline(y = 0, color = color, lw = 5)
        ax[subsection[i]].set_ylabel("")
        ax[subsection[i]].sharex(ax["B"])
        ax[subsection[i]].set_xticks(ticks=[],
                           labels=[])

    if subsection[i] == subsection[-1]:
        ax[subsection[i]].set_xlabel("Estimated direction of age-dependency", fontsize=9)
    # axes[i, j].set_titles("")
    ax[subsection[i]].set(yticks=[])
    # axes[i, j].set_ylabels()
    # axes[i, j].despine(left=True, bottom=True)
    ax[subsection[i]].axvline(x=1, linestyle="dashed", color="black")

    sns.despine(left=True, bottom= True)
    plt.tight_layout()