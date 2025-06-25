"""
Project: Diversification Rates and ADE
Description:
Extraction of Weibull shape parameters predicted by ADE-Bayes, calculation of mean + 95% CIs
"""

from pandas import *
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns

path = "/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2024/ADE-Bayes/new_PyRate"

time_ints= [-95.61, -83.5, -73.1, -71.99, -66.59, -65.79, -55.76, -38.77, -33.47, -3.45]

weibul = []
ci_min = []
ci_max = []

for i in range(1, 11):
    data = read_csv(path + "/all_species_{i}_ADE_ADE_mcmc.log".format(i  = i), sep = "\t")

    weibul.append(np.mean(data["w_shape"]))
    ste = st.sem(data["w_shape"])
    ci_min.append(data['w_shape'].quantile(0.025))
    ci_max.append(data['w_shape'].quantile(0.975))

# plot

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

sns.scatterplot(x = time_ints, y = weibul, ax = ax["A"])
ax["A"].hlines(y = 1, xmin=-145, xmax = 0, color = "black", linestyles= "dashed")
ax["A"].vlines(x = time_ints, ymin = ci_min, ymax = ci_max, color = "black")

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
