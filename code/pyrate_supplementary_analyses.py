"""
Project: Diversification Rates and ADE
Description:
Additional visualisation of raw data in combination with PyRate estimated rates
1. Number of min_mas and extinction rates + Number of max_mas and origination rates (check whether accummulation of these occurr at the same time as high rates)
2. Number of collections per stage (A "simple" way of looking at preservation)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import matplotlib.patheffects as pe
from collections import Counter

def convert(object):
    lyst = object.to_list()
    for i in lyst:
        lyst = i.split(",")
    lyst[-1] = lyst[-1].replace(")", "")
    lyst[0] = re.sub(".*" + "\(", "", lyst[0])
    for i in range(len(lyst)):
        lyst[i] = float(lyst[i])
    return lyst

occs = pd.read_csv("/Users/kristinakocakova/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/inputs/2025/species.txt", sep = "\t")
occs["max_ma"] = occs["max_ma"]*-1
occs["min_ma"] = occs["min_ma"]*-1

count_min = dict(sorted(Counter(occs["min_ma"]).items()))
count_max = dict(sorted(Counter(occs["max_ma"]).items()))

min = list(count_min.keys())
max = list(count_max.keys())

no_min = list(count_min.values())
no_max = list(count_max.values())


rates = pd.read_csv("/Users/kristinakocakova/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species/RTT_plots.r", sep="\t", header=None)

time_e = convert(rates.iloc[19])
rate_e = convert(rates.iloc[20])


time_s = convert(rates.iloc[3])
rate_s = convert(rates.iloc[4])

# plot Min_ma and ext rate

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
sns.set_style("whitegrid")

fig, ax = plt.subplot_mosaic(
    """
    AA
    BB
    CC
    """,
    figsize=(6.69291, 4), gridspec_kw={"height_ratios": [8, 1, 0.5]}, sharex = True
)

ax["A"].bar(min, [x/1000 for x in no_min], alpha=0.3, color="purple", zorder=1)
sns.lineplot(x = time_e, y = rate_e, color="#bf3e36", ax = ax["A"], linewidth=1.5,
             path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()], zorder = 2)

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

ax["B"].bar(stage_chrono, bar_vals_stage, width=width_stage, color=cols_stage, edgecolor="none")
ax["B"].bar_label(ax["B"].containers[0], labels=stage_chrono_names, padding=-15, fontsize=6, rotation=90)
ax["B"].set_xlim(xmin=-146)
ax["B"].set_ylim(ymin=0, ymax=0.4)
ax["B"].set_yticks(ticks=[0], labels=[""])
ax["B"].axis("off")

# Time-line - Epochs

epoch_chrono = [-122.5, -83, -61, -44.95, -28.465, -14.1815, -3.9565, -1.29585, -0.00585]
bar_vals = [0.25] * 9
width = [44.5, 34, 10, 22.1, 10.9, 17.7, 2.753, 2.5683, 0.0117]
cols = ["#8fd07b", "#67a155", "#ffb17b", "#ffbc87", "#ffc694", "#ffeb3e", "#fff6b2", "#ffefaf", "#feebd2"]
epoch_chrono_names = ["Early Cretaceous", "Late Cretaceous", "Paleoc.", "Eocene", "Oligocene", "Miocene", "Pli.",
                      "Ple.", ""]

ax["C"].bar(epoch_chrono, bar_vals, width=width, color=cols, edgecolor="none")
ax["C"].bar_label(ax["C"].containers[0], labels=epoch_chrono_names, padding=-10, fontsize=6)
ax["C"].set_xlim(xmin=-146)
ax["C"].set_ylim(ymin=0, ymax=0.25)
ax["C"].set_yticks(ticks=[0], labels=[""])

ax["C"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])
ax["C"].set_xlabel("Time (Ma)", fontsize=7)
ax["C"].axis("off")
ax["C"].set_visible(True)

plt.tight_layout()
sns.despine()



# plot Max_ma and orig rate

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
sns.set_style("whitegrid")

fig, ax = plt.subplot_mosaic(
    """
    AA
    BB
    CC
    """,
    figsize=(6.69291, 4), gridspec_kw={"height_ratios": [8, 1, 0.5]}, sharex = True
)

ax["A"].bar(max, [x/1000 for x in no_max], alpha=0.3, color="purple", zorder=1)
sns.lineplot(x = time_s, y = rate_s, color="#2d7096", ax = ax["A"], linewidth=1.5,
             path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()], zorder = 2)

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

ax["B"].bar(stage_chrono, bar_vals_stage, width=width_stage, color=cols_stage, edgecolor="none")
ax["B"].bar_label(ax["B"].containers[0], labels=stage_chrono_names, padding=-15, fontsize=6, rotation=90)
ax["B"].set_xlim(xmin=-146)
ax["B"].set_ylim(ymin=0, ymax=0.4)
ax["B"].set_yticks(ticks=[0], labels=[""])
ax["B"].axis("off")

# Time-line - Epochs

epoch_chrono = [-122.5, -83, -61, -44.95, -28.465, -14.1815, -3.9565, -1.29585, -0.00585]
bar_vals = [0.25] * 9
width = [44.5, 34, 10, 22.1, 10.9, 17.7, 2.753, 2.5683, 0.0117]
cols = ["#8fd07b", "#67a155", "#ffb17b", "#ffbc87", "#ffc694", "#ffeb3e", "#fff6b2", "#ffefaf", "#feebd2"]
epoch_chrono_names = ["Early Cretaceous", "Late Cretaceous", "Paleoc.", "Eocene", "Oligocene", "Miocene", "Pli.",
                      "Ple.", ""]

ax["C"].bar(epoch_chrono, bar_vals, width=width, color=cols, edgecolor="none")
ax["C"].bar_label(ax["C"].containers[0], labels=epoch_chrono_names, padding=-10, fontsize=6)
ax["C"].set_xlim(xmin=-146)
ax["C"].set_ylim(ymin=0, ymax=0.25)
ax["C"].set_yticks(ticks=[0], labels=[""])

ax["C"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])
ax["C"].set_xlabel("Time (Ma)", fontsize=7)
ax["C"].axis("off")
ax["C"].set_visible(True)

plt.tight_layout()
sns.despine()