"""
Project: Diversification Rates and ADE
Description:
Additional visualisation of raw data in combination with PyRate estimated rates to assess potential biases in PyRate
Based on recommendations from https://academic.oup.com/sysbio/article/71/1/153/6295892?login=true
1. Number of min_mas and extinction rates + Number of max_mas and origination rates (check whether accummulation of these occurr at the same time as high rates)
2. Number of collections per stage (A "simple" way of looking at preservation)
3. Number of occurrences per stage
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import matplotlib.patheffects as pe
from collections import Counter

# def convert(object):
#     lyst = object.to_list()
#     for i in lyst:
#         lyst = i.split(",")
#     lyst[-1] = lyst[-1].replace(")", "")
#     lyst[0] = re.sub(".*" + "\(", "", lyst[0])
#     for i in range(len(lyst)):
#         lyst[i] = float(lyst[i])
#     return lyst

occs = pd.read_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K_version/Manuscripts/Rates_+_ADE/For submission/Proceedings B/R1/code and data/Inputs/PyRate/species.txt", sep = "\t")

# occs["max_ma"] = occs["max_ma"]*-1
# occs["min_ma"] = occs["min_ma"]*-1
#
# count_min = dict(sorted(Counter(occs["min_ma"]).items()))
# count_max = dict(sorted(Counter(occs["max_ma"]).items()))
#
# min = list(count_min.keys())
# max = list(count_max.keys())
#
# no_min = list(count_min.values())
# no_max = list(count_max.values())
#
#
# rates = pd.read_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species/RTT_plots.r", sep="\t", header=None)
#
# time_e = convert(rates.iloc[19])
# rate_e = convert(rates.iloc[20])
#
#
# time_s = convert(rates.iloc[3])
# rate_s = convert(rates.iloc[4])

# 2. calculate number of collections per stage

# modify the dataset to contain only one representative of a collection per time bin

df_cols = occs.drop_duplicates(subset = ["max_ma", "min_ma", "collection_no"])

time_ints = [145, 139.8, 132.6, 125.77, 121.4, 113, 100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66, 61.6, 59.2, 56, 47.8, 41.2, 37.71, 33.9, 27.82, 23.03, 20.44, 15.98, 13.82, 11.63, 7.246,
                5.333, 3.6, 2.58, 1.8, 0.774, 0.129, 0.0117, 0]

cols_count = []

for i in range(len(time_ints) - 1):
    NFL = len(df_cols.loc[(df_cols["max_ma"] <= time_ints[i]) & (df_cols["min_ma"] >= time_ints[i + 1])])
    NbL = len(df_cols.loc[(df_cols["max_ma"] > time_ints[i]) & (df_cols["min_ma"] > time_ints[i + 1]) & (df_cols["min_ma"] < time_ints[i])])
    NFt = len(df_cols.loc[(df_cols["max_ma"] <= time_ints[i]) & (df_cols["max_ma"] > time_ints[i + 1]) & (df_cols["min_ma"] < time_ints[i + 1])])
    Nbt = len(df_cols.loc[(df_cols["max_ma"] > time_ints[i]) & (df_cols["min_ma"] < time_ints[i + 1])])
    Ntot = NFL + NbL + NFt + Nbt
    cols_count.append(Ntot)

# plot Min_ma and ext rate
time_ints = [x* -1 for x in time_ints]

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
sns.set_style("white")

fig, ax = plt.subplot_mosaic(
    """
    BB
    CC
    DD
    EE
    """,
    figsize=(6.69291, 6), gridspec_kw={"height_ratios": [4, 4, 1, 0.5]}, sharex = True
)

fig.subplots_adjust(hspace=0.05)

ax["C"].step(x = time_ints[:-1], y = cols_count)
ax["C"].set_ylabel("Number of collections\nper stage", fontsize=7)


# Time-line - Stages

stage_chrono = [-142.4, -136.2, -129.185, -123.59, -117.2, -106.75, -97.2, -91.85, -88.05, -84.95, -77.85, -69.05,
                -63.8, -60.4, -57.6, -51.9, -44.5, -39.455, -35.805, -30.86, -25.425, -21.735, -18.205,
                -14.895, -12.725, -9.438, -6.2895, -3.9565, -1.29585, -0.00585]
bar_vals_stage = [0.4] * 30
width_stage = [5.2, 7.2, 6.83, 4.37, 8.4, 12.5, 6.6, 4.1, 3.5, 2.7, 11.5, 6.1, 4.4, 2.4, 3.2, 8.2, 6.6, 3.5, 3.81, 6.08,
               4.79, 2.59, 4.47, 2.15, 2.19, 4.384, 1.913,
               2.753, 2.5683, 0.0117]
cols_stage = ["#8ed084", "#9ad58d", "#a7d996", "#b4de9f", "#bfe48a", "#c1e2a8", "#cee7b1", "#b8da7a", "#c5df82",
              "#d2e38c", "#dfe895", "#eaec9e", "#fdb462", "#ffc482",
              "#ffc58b", "#ffb38b", "#ffbe98", "#ffc8a5", "#ffd2b3", "#ffdbae", "#ffe5bc", "#ffee65", "#ffef6e",
              "#fff078", "#fff181", "#fff28b", "#fff395", "#fff6b2",
              "#ffefaf", "#feebd2"]
stage_chrono_names = ["Ber.", "Val.", "Hau.", "Bar.", "Apt.", "Alb.", "Cen.", "Tur.", "Con.", "San.", "Cam.",
                      "Maa.", "Dan.", "Sel.", "Tha.", "Ypr.", "Lut.", "Bart.", "Pri.", "Rup.", "Cha.", "Aqu.",
                      "Bur.", "Lan.", "Ser.", "Tor.", "Mes.", "", "", ""]

ax["D"].bar(stage_chrono, bar_vals_stage, width=width_stage, color=cols_stage, edgecolor="none")
ax["D"].bar_label(ax["D"].containers[0], labels=stage_chrono_names, padding=-15, fontsize=6, rotation=90)
ax["D"].set_xlim(xmin=-146)
ax["D"].set_ylim(ymin=0, ymax=0.4)
ax["D"].set_yticks(ticks=[0], labels=[""])
ax["D"].axis("off")

# Time-line - Epochs

epoch_chrono = [-122.5, -83, -61, -44.95, -28.465, -14.1815, -3.9565, -1.29585, -0.00585]
bar_vals = [0.25] * 9
width = [44.5, 34, 10, 22.1, 10.9, 17.7, 2.753, 2.5683, 0.0117]
cols = ["#8fd07b", "#67a155", "#ffb17b", "#ffbc87", "#ffc694", "#ffeb3e", "#fff6b2", "#ffefaf", "#feebd2"]
epoch_chrono_names = ["Early Cretaceous", "Late Cretaceous", "Paleoc.", "Eocene", "Oligocene", "Miocene", "Pli.",
                      "Ple.", ""]

ax["E"].bar(epoch_chrono, bar_vals, width=width, color=cols, edgecolor="none")
ax["E"].bar_label(ax["E"].containers[0], labels=epoch_chrono_names, padding=-10, fontsize=6)
ax["E"].set_xlim(xmin=-146)
ax["E"].set_ylim(ymin=0, ymax=0.25)
ax["E"].set_yticks(ticks=[0], labels=[""])

ax["E"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])
ax["E"].set_xlabel("Time (Ma)", fontsize=7)
ax["E"].axis("off")
ax["E"].set_visible(True)

sns.despine()
plt.tight_layout()


# 3. plot of occurrence number per epoch

data = pd.ExcelFile("/Volumes/External_memory/Dropbox/Kristina_PhD_K_version/Manuscripts/Rates_+_ADE/For submission/Proceedings B/R1/code and data/species+genera.xlsx")
occs = pd.read_excel(data, sheet_name= "Species") # or Genera
# 329 occurrences have different early and late epoch

epochs = list(occs["early_epoch"].value_counts().keys())
epoch_counts = list(occs["early_epoch"].value_counts())

# epoch_counts_scaled = []
#
# epoch_durs = {"Lower Cretaceous": 44.5, "Upper Cretaceous": 34.5, "Paleocene": 10, "Eocene": 22.1, "Oligocene": 10.87, "Miocene": 17.697, "Pliocene": 2.753, "Pleistocene": 2.5683}
#
# for i, j in enumerate(epochs):
#     epoch_counts_scaled.append(epoch_counts[i]/epoch_durs[j])

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

fig, ax = plt.subplots(figsize = (6.69291, 8.26772/2))

sns.barplot(x = epochs, y = epoch_counts,  hue = epochs, palette= {"Lower Cretaceous": "#8CCD57", "Upper Cretaceous": "#A6D84A", "Paleocene": "#FDA75F", "Eocene": "#FDB46C", "Oligocene": "#FDC07A", "Miocene": "#FFFF00", "Pliocene": "#FFFF99", "Pleistocene": "#FFF2AE", "Holocene" : "#ffffff"},
                order= ["Lower Cretaceous", "Upper Cretaceous","Paleocene", "Eocene", "Oligocene", "Miocene", "Pliocene","Pleistocene", "Holocene"])

for container in ax.containers:
    ax.bar_label(container, color = "#3a3a3a")
ax.set_xticklabels(["Lower Cret.", "Upper Cret.","Paleocene", "Eocene", "Oligocene", "Miocene", "Pliocene","Pleistocene", "Holocene"])
plt.ylabel("Number of occurrences")
plt.xlabel("Epoch")
sns.despine()

plt.show(block = True)





