"""
Project: Diversification Rates and ADE
Description:
Plotting of a skyline plot of the patterns obtained using the DivDyn package
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

d = "/Users/kristinakocakova/" # or /Volumes/External_memory/

df = pd.read_csv(d + "Dropbox/Kristina_PhD_K_version/Kristina_files/Analyses/PyRate/PyRate_Analysis/inputs/2025/June/div_dyn_rates.txt", sep = "\t")

df['max_ma'] *= -1

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
    DD
    EE
    """,
    figsize=(6.69291, 6), gridspec_kw={"height_ratios": [3, 3, 3, 1, 0.5]}, sharex = True
)

fig.subplots_adjust(hspace=0.05)

ax["A"].step(x = df["max_ma"], y = df["extPC_mean"], color = "#bf3e36")
ax["A"].fill_between(df["max_ma"], df["extPC_lower"], df["extPC_upper"], step="pre", color='#bf3e36', alpha=0.15)
ax["A"].set_ylabel("Per capita\nextinction rate", fontsize=7)

ax["B"].step(x = df["max_ma"], y = df["oriPC_mean"], color = "#2d7096")
ax["B"].fill_between(df["max_ma"], df["oriPC_lower"], df["oriPC_upper"], step="pre", color='#2d7096', alpha=0.15)
ax["B"].set_ylabel("Per capita\nspeciation rate", fontsize=7)

ax["C"].step(x = df["max_ma"], y = df["net_div_mean"], color = "#4f5450")
ax["C"].fill_between(df["max_ma"], df["net_div_lower"], df["net_div_upper"], step="pre", color='#2d7096', alpha=0.15)
ax["C"].set_ylabel("Per capita\nnet diversification rate", fontsize=7)


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



