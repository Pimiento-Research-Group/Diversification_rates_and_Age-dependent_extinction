"""
Project: Diversification Rates and ADE
Description:
Calculate and plot extinction, speciation and net diversification based on FADs and LADs
Total diversity (Ntot) = NFL + NFt + NbL + Nbt
Number of originations = NFL + NFt
Number of extinctions = NFL + NbL
Origination rate = (NFL + NFt)/Ntot/t
Extinction rate = (NFL + NbL)/Ntot/t
Net diversification rate = origination - extinction

NFL - confined to the interval - FAD < t0 & LAD > t1
NbL - bottom boundary crossers - FAD > t0 & LAD > t1 & LAD < t0
NFt - top boundary crossers - FAD < t0 & FAD > t1 & LAD < t1
Nbt - both boundaries crossed - FAD > t0 & LAD < t1

https://www.jstor.org/stable/1571654?seq=3

Plot preservation rate through time (q_rate from combined_mcmc.log)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patheffects as pe


species = pd.read_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species_raw.txt", sep = "\t")

time_ints = [145, 139.8, 132.6, 125.77, 121.4, 113, 100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66, 61.6, 59.2, 56, 47.8, 41.2, 37.71, 33.9, 27.82, 23.03, 20.44, 15.98, 13.82, 11.63, 7.246,
                5.333, 3.6, 2.58, 1.8, 0.774, 0.129, 0.0117, 0]

ext_rate = []
orig_rate = []
net_rate = []
ext_no = []
orig_no = []
n_tot = []

for i in range(len(time_ints) - 1):
    NFL = len(species.loc[(species["fad"] <= time_ints[i]) & (species["lad"] >= time_ints[i + 1])])
    NbL = len(species.loc[(species["fad"] >= time_ints[i]) & (species["lad"] >= time_ints[i + 1]) & (species["lad"] <= time_ints[i])])
    NFt = len(species.loc[(species["fad"] <= time_ints[i]) & (species["fad"] >= time_ints[i + 1]) & (species["lad"] <= time_ints[i + 1])])
    Nbt = len(species.loc[(species["fad"] >= time_ints[i]) & (species["lad"] <= time_ints[i + 1])])
    Ntot = NFL + NbL + NFt + Nbt

    ext_no.append(NFL + NbL)
    orig_no.append(NFL + NFt)
    e = (NFL + NbL)/Ntot
    o = (NFL + NFt)/Ntot
    ext_rate.append(e)
    orig_rate.append(o)
    net_rate.append(o - e)
    n_tot.append(Ntot)

df_rates = pd.DataFrame()
df_rates["stage"] = time_ints[:-1]
df_rates["ext_rate"] = ext_rate
df_rates["orig_rate"] = orig_rate
df_rates["net"] = net_rate
df_rates["ext_no"] = ext_no
df_rates["orig_no"] = orig_no
df_rates["n_tot"] = n_tot


# plot the rates

# convert time ints to match
df_rates["stage"] = df_rates["stage"]*-1

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

fig, ax = plt.subplot_mosaic(
    """
    AA
    BB
    CC
    DD
    EE
    """,
    figsize=(6.69291, 6), gridspec_kw={"height_ratios": [3, 3, 3, 1, 0.5]}, sharex = True
)

sns.despine()
plt.tight_layout()
fig.subplots_adjust(hspace=0.05)
ax["E"].xaxis.set_tick_params(which="both", labelbottom=True)
sns.set_style("whitegrid")

# Extinction
sns.lineplot(data = df_rates, x = "stage", y = "ext_rate", color="#bf3e36", ax = ax["A"], linewidth=1.5,
             path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()], zorder = 1)
sns.scatterplot(data = df_rates, x = "stage", y = "ext_rate", ax = ax["A"], color="#bf3e36", zorder = 2)
ax["A"].set_ylim(ymin=-0.3, ymax=0.7)
ax["A"].set_xlim(xmin = -146, xmax = 1)
ax["A"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])
ax["A"].set_ylabel("Extinction rate", fontsize=7)

# Speciation

sns.lineplot(data = df_rates, x = "stage", y = "orig_rate", ax = ax["B"], color="#2d7096", linewidth=1.5,
             path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()], zorder = 1)
sns.scatterplot(data = df_rates, x = "stage", y = "orig_rate", ax = ax["B"] , color="#2d7096", zorder = 2)
ax["B"].set_ylim(ymin=-0.3, ymax=1.1)
ax["B"].set_xlim(xmin = -146, xmax = 1)
ax["B"].set_yticks(ticks=[0, 0.5, 1.0],
                   labels=[0, 0.5, 1.0])
ax["B"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])
ax["B"].set_ylabel("Speciation rate", fontsize=7)

# Net Diversification

sns.lineplot(data = df_rates, x = "stage", y = "net", ax = ax["C"], color="#4f5450", linewidth=1.5,
             path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()], zorder = 1)
sns.scatterplot(data = df_rates, x = "stage", y = "net", ax = ax["C"], color="#4f5450", zorder = 2)
ax["C"].set_ylim(ymin=-1, ymax=1.1)
ax["C"].set_xlim(xmin = -146, xmax = 1)
ax["C"].set_ylabel("Net diversification rate", fontsize=7)
ax["C"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])

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


plt.tight_layout()
plt.show()







