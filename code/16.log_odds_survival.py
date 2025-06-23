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
import statsmodels.api as sm
import statsmodels.formula.api as smf

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
probs = []
bin = []

for i in time_int_max:
    tb = df_bins.loc[df_bins["bin"] == i]
    print(len(tb))
    if len(tb) > 1 and len(np.unique(tb["survival"].values)) > 1:
        lgit =smf.logit("survival ~ age", data = tb)
        res= lgit.fit(disp = False)
        lgits.append(res.params["age"])
        probs.append(res.pvalues["age"])
        bin.append(i)


df_odds = pd.DataFrame()
df_odds["bin"] = bin
df_odds["log_reg"] = lgits
df_odds["prob"] = probs
df_odds["mid-bin"] = df_odds["bin"] - 2.5
df_odds["mid-bin"] = df_odds["mid-bin"] * -1

sig = []
for i in df_odds["prob"].values:
    if i < 0.05:
        sig.append("yes")
    else:
        sig.append("no")

df_odds["sig"] = sig

# plot these through time

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

fig, ax = plt.subplot_mosaic(
    """
    AA
    BB
    CC
    """,
    figsize=(6.69291, 4), gridspec_kw={"height_ratios": [3, 0.5, 0.5]}, sharex=True
)

sns.despine()
plt.tight_layout()
fig.subplots_adjust(hspace=0.05)
ax["B"].xaxis.set_tick_params(which="both", labelbottom=True)

sns.scatterplot(data = df_odds.loc[df_odds["bin"] != 55], x = "mid-bin", y = "log_reg", hue = "sig", ax = ax["A"])
ax["A"].hlines(y = 0, xmin=-145, xmax = 0, color = "black", linestyles= "dashed")

# Timeline - Stages

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
ax["B"].set_xlim(xmin=-145)
ax["B"].set_ylim(ymin=0, ymax=0.4)
ax["B"].set_yticks(ticks=[0], labels=[""])
ax["B"].axis("off")
# ax["D"].set_xticklabels([])

# Timeline - Epochs

epoch_chrono = [-122.5, -83, -61, -44.95, -28.465, -14.1815, -3.9565, -1.29585, -0.00585]
bar_vals = [0.25] * 9
width = [44.5, 34, 10, 22.1, 10.9, 17.7, 2.753, 2.5683, 0.0117]
cols = ["#8fd07b", "#67a155", "#ffb17b", "#ffbc87", "#ffc694", "#ffeb3e", "#fff6b2", "#ffefaf", "#feebd2"]
epoch_chrono_names = ["Early Cretaceous", "Late Cretaceous", "Paleoc.", "Eocene", "Oligocene", "Miocene", "Pli.",
                      "Ple.", ""]

ax["C"].bar(epoch_chrono, bar_vals, width=width, color=cols, edgecolor="none")
ax["C"].bar_label(ax["C"].containers[0], labels=epoch_chrono_names, padding=-10, fontsize=6)
ax["C"].set_xlim(xmin=-145)
ax["C"].set_ylim(ymin=0, ymax=0.25)
ax["C"].set_yticks(ticks=[0], labels=[""])

ax["C"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])
ax["C"].set_xlabel("Time (Ma)", fontsize=7)
ax["C"].axis("off")
ax["C"].set_visible(True)

plt.show()





