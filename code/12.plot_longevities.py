"""
Project: Diversification Rates and ADE
Description:
Plotting of figure S7 - longevities distribution of extinct taxa
Input file - .txt file generated using the -ginput function inscript 3.2.
"""

from pandas import *
options.display.width = 0
import matplotlib.pyplot as plt
from statistics import mean, median
import seaborn as sns
import numpy as np


#names of files
lst = ["species", "genera", "species_raw", "genera_raw"]
#axes
indx_lst = ["A", "B", "C", "D"]

p = "/path_to_pyrate_output_folder/"

#/Users/kristinakocakova/Dropbox/Kristina_PhD_K_version/Kristina_files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6
# fig, ax = plt.subplots(2 , 2, figsize = (6.69291, 8.26772/4))
fig, ax = plt.subplot_mosaic(
    """
    AB
    CD
    """, figsize = (6.69291, 8.26772/2)
)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)

for i in range(len(lst)):
    if "raw" in lst[i]:
        te_ts_sp = read_csv(p + "{i}.txt".format(i=lst[i]), sep="\t")

        te_ts_sp = te_ts_sp.loc[te_ts_sp["lad"] != 0]

        # for time bins
        # te_ts_sp = te_ts_sp.loc[(te_ts_sp["te"] <= 66.78)&(te_ts_sp["te"] >= 65.38)]

        te = te_ts_sp["lad"].to_numpy()
        ts = te_ts_sp["fad"].to_numpy()

        longevity = ts - te

        longevity = longevity.tolist()

        mean_long = np.mean(longevity)

    else:
        te_ts_sp = read_csv(p + "{i}_se_est.txt".format(i=lst[i]), sep="\t")

        te_ts_sp = te_ts_sp.loc[te_ts_sp["te"] != 0]

        # for time bins
        # te_ts_sp = te_ts_sp.loc[(te_ts_sp["te"] <= 66.78)&(te_ts_sp["te"] >= 65.38)]

        te = te_ts_sp["te"].to_numpy()
        ts = te_ts_sp["ts"].to_numpy()

        longevity = ts - te

        longevity = longevity.tolist()


        mean_long = np.mean(longevity)

    sns.set_style("whitegrid")
    h = ax[indx_lst[i]].hist(longevity, edgecolor = "white", color = "#a8293c", bins = 20)
    ymax = np.max(h[0])
    sns.despine()
    # sns.despine(ax=ax[0, indx_lst[i]])
    # sns.despine(ax=ax[0, indx_lst[i]])
    ax[indx_lst[i]].vlines(x=[mean_long], ymin=0, ymax = ymax, linestyles="dashed", colors=['#000000'], alpha=0.5, linewidth = 2)
    ax[indx_lst[i]].text(mean_long + 3, ymax - 5, s =str(round(mean_long, 3)), c = '#000000', fontsize = 9)
    # plt.yticks(fontsize=15)
    # plt.xticks(fontsize=15)
    # ax1.tick_params(axis='both', which='major', labelsize=7)

plt.subplots_adjust(hspace=0.5)
fig.supxlabel("Longevity (Myr)", fontsize=10)
fig.supylabel("Frequency", fontsize=10)
plt.tight_layout()

