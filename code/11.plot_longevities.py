"""
Project: SHARK-XT Rates and ADE
Author: Kristína Kocáková
Description:
Plotting of figure S7 - longevities distribution
"""

from pandas import *
options.display.width = 0
import matplotlib.pyplot as plt
from statistics import mean, median
import seaborn as sns
import numpy as np


#names of files
lst = ["species", "selachi", "batos", "lamni", "carcharhini", "orectolobi", "hexanchi", "squali", "mylio", "raji", "rhinopristi"]
#axes
indx_lst = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1]

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6
fig, ax = plt.subplots(4 , 3, figsize = (6.69291, 8.26772), sharex=True)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)

for i in range(len(lst)):
    if i < 3:
        ax1 = ax[0, indx_lst[i]]
    elif i >=3 and i <6:
        ax1 = ax[1, indx_lst[i]]
    elif i >=6 and i <9:
        ax1 = ax[2, indx_lst[i]]
    elif i >=9 and i <12:
        ax1 = ax[3, indx_lst[i]]
    else:
        ax1 = ax[4, indx_lst[i]]
    te_ts_sp = read_csv("/Users/kristinakocakova/Dropbox/Kristina_PhD/Analyses/PyRate/PyRate_Analysis/outputs/2024/fast_burnin_10_2024/{i}/se_est.txt".format(i=lst[i]), sep="\t")

    te_ts_sp = te_ts_sp.loc[te_ts_sp["te"] != 0]

    # for time bins
    # te_ts_sp = te_ts_sp.loc[(te_ts_sp["te"] <= 66.78)&(te_ts_sp["te"] >= 65.38)]

    te = te_ts_sp["te"].to_numpy()
    ts = te_ts_sp["ts"].to_numpy()

    longevity = ts - te

    longevity = longevity.tolist()


    mean_long = np.mean(longevity)

    sns.set_style("whitegrid")
    h = ax1.hist(longevity, edgecolor = "white", color = "#a8293c", bins = 20)
    ymax = np.max(h[0])
    sns.despine()
    # sns.despine(ax=ax[0, indx_lst[i]])
    # sns.despine(ax=ax[0, indx_lst[i]])
    ax1.vlines(x=[mean_long], ymin=0, ymax = ymax, linestyles="dashed", colors=['#000000'], alpha=0.5, linewidth = 2)
    ax1.text(mean_long + 3, ymax - 5, s =str(round(mean_long, 3)), c = '#000000', fontsize = 9)
    # plt.yticks(fontsize=15)
    # plt.xticks(fontsize=15)
    # ax1.tick_params(axis='both', which='major', labelsize=7)
    plt.tight_layout()
    plt.show()

fig.supxlabel("Longevity (Myr)", fontsize=10)
fig.supylabel("Frequency", fontsize=10)

