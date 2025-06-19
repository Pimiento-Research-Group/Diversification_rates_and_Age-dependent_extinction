"""
Project: Diversification Rates and ADE
Description:
Function to produce a plot of extinction, speciation and net diversification, and to save the estimated rates in an .xlsx file
Input file - .R script generated using the -plotRJ
Figure 1. - Plotted without rate shifts, simply remove or comment out FG from the plt.subplot.mosaic and the ax lines with ["F"] or ["G"]
Figure S2-S3 - Plotted with the full plot() function
Figure S1 - Individual rate shifts, see below plot() function
"""

from pandas import *
import re
import matplotlib.pyplot as plt
import scipy.stats as st
import numpy as np
import seaborn as sns
import matplotlib.patheffects as pe
plt.rcParams['svg.fonttype'] = 'none'

def plot(file_name):

    r = read_csv("/path_to_pyrate_output_folder/{i}/RTT_plots.r".format(i = file_name), sep="\t", header= None)

    # CONVERT PYRATE GENERATED R FILE INTO PYTHON LISTS OF SPECIATION AND EXTINCTION RATE VALUES
    def convert(object):
        lyst = object.to_list()
        for i in lyst:
            lyst = i.split(",")
        lyst[-1] = lyst[-1].replace(")", "")
        lyst[0] = re.sub(".*" + "\(", "", lyst[0])
        for i in range(len(lyst)):
            lyst[i] = float(lyst[i])
        return lyst
    def convert_bf(object):
        lyst = object.to_list()
        for i in lyst:
            x = i.split()
        y = float(x[2])
        return y
    """
    IMPORTANT
    The row numbers may differ depending on the output file from PyRate, which in turn depends on the parameters.
    So it's best to first check in the .r file where exactly rates, CIs and times are stored. 
    """
    #EXTINCTION
    time_e = convert(r.iloc[19])
    rate_e = convert(r.iloc[20])
    minHPD_e = convert(r.iloc[21])
    maxHPD_e = convert(r.iloc[22])

    #SPECIATION
    time_s = convert(r.iloc[3])
    rate_s = convert(r.iloc[4])
    minHPD_s = convert(r.iloc[5])
    maxHPD_s = convert(r.iloc[6])

    #NET DIVERSIFICATION
    time_d = convert(r.iloc[34])
    rate_d = convert(r.iloc[35])
    minHPD_d = convert(r.iloc[36])
    maxHPD_d = convert(r.iloc[37])

    #Extinction rate shifts
    time_es = convert(r.iloc[28])
    freq_es = convert(r.iloc[29])
    bf2 = convert_bf(r.iloc[15])
    bf6 = convert_bf(r.iloc[16])

    #Speciation rate shifts
    time_ss = convert(r.iloc[12])
    freq_ss = convert(r.iloc[13])

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
        FG
        """,
        figsize=(6.69291, 6), gridspec_kw={"height_ratios":[3, 3, 3, 1, 0.5, 4]}, sharex=True
    )

    sns.despine()
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    ax["E"].xaxis.set_tick_params(which= "both",labelbottom=True)

    # max_y_e = max(maxHPD_e)
    sns.set_style("whitegrid")
    #Extinction
    ax["A"].plot(time_e, rate_e, color="#bf3e36", linewidth = 1.5, path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()])
    ax["A"].fill_between(
        time_e, minHPD_e, maxHPD_e, color='#bf3e36', alpha=.15)
    ax["A"].set_ylim(ymin=0, ymax = 1.3)

    # harm_mean_e = len(rate_e) / np.sum(1 / np.array(rate_e))
    # ax["A"].plot([-145, 0], [harm_mean_e, harm_mean_e], color="#038249",  linewidth = 1, linestyle = "dashed", path_effects=[pe.Stroke(linewidth=3, foreground='w', alpha = .7), pe.Normal()], alpha =.7)
    ax["A"].set_ylim(ymin=0, ymax = 1.3)
    ax["A"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
               labels=["", "", "", "", "", ""])

    # ax["A"].set_xticklabels([])

    ax["A"].set_ylabel("Extinction rate", fontsize = 7)
    legend_elements = [plt.Line2D([0], [0], color='#038249', lw=1, label='Background rate', linestyle= "dashed")]
    ax["A"].legend(handles = legend_elements, loc = "upper right", fontsize = 7)

    #Speciation

    ax["B"].plot(time_s, rate_s, color="#2d7096", linewidth = 1.5, path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()])
    ax["B"].fill_between(
        time_s, minHPD_s, maxHPD_s, color='#2d7096', alpha=.15)
    ax["B"].set_ylim(ymin=0, ymax = 1.1)

    # harm_mean_s = len(rate_e) / np.sum(1 / np.array(rate_e))
    # ax["B"].plot([-145, 0], [harm_mean_s, harm_mean_s], color="#038249",  linewidth = 1.5, linestyle = "dashed", path_effects=[pe.Stroke(linewidth=3.5, foreground='w', alpha = .7), pe.Normal()], alpha =.7)
    ax["B"].set_ylim(ymin=0, ymax = 1.1)
    ax["B"].set_yticks(ticks=[0, 0.5, 1.0],
               labels=[0, 0.5, 1.0])
    ax["B"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
               labels=["", "", "", "", "", ""])


    # ax["B"].set_xticklabels([])

    ax["B"].set_ylabel( "Speciation rate", fontsize = 7)

    # Net Diversification

    ax["C"].plot(time_d, rate_d, color = "#4f5450", linewidth = 1.5, path_effects=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()])
    ax["C"].fill_between(
        time_d, minHPD_d, maxHPD_d, color='grey', alpha=.4)
    ax["C"].set_ylim(ymin=-1, ymax = 1.1)

    # ax["C"].plot([-143, 0], [0, 0], linestyle="dashed", color = "black")

    ax["C"].set_ylabel("Net diversification rate", fontsize = 7)
    ax["C"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
               labels=["", "", "", "", "", ""])

    # ax["C"].set_xticklabels([])

    #Time line - Stages

    stage_chrono = [-142.4, -136.2, -129.185, -123.59, -117.2, -106.75, -97.2, -91.85, -88.05, -84.95, -77.85, -69.05, -63.8, -60.4, -57.6, -51.9, -44.5, -39.455, -35.805, -30.86, -25.425, -21.735, -18.205,
                    -14.895, -12.725, -9.438, -6.2895, -3.9565, -1.29585, -0.00585]
    bar_vals_stage =  [0.4]*30
    width_stage = [5.2, 7.2, 6.83, 4.37, 8.4, 12.5, 6.6, 4.1, 3.5, 2.7, 11.5, 6.1, 4.4, 2.4, 3.2, 8.2, 6.6, 3.5, 3.81, 6.08, 4.79, 2.59, 4.47, 2.15, 2.19, 4.384, 1.913,
                   2.753, 2.5683, 0.0117]
    cols_stage = ["#8ed084", "#9ad58d", "#a7d996", "#b4de9f", "#bfe48a", "#c1e2a8", "#cee7b1", "#b8da7a", "#c5df82", "#d2e38c", "#dfe895", "#eaec9e", "#fdb462", "#ffc482",
                  "#ffc58b", "#ffb38b", "#ffbe98", "#ffc8a5", "#ffd2b3", "#ffdbae", "#ffe5bc", "#ffee65", "#ffef6e", "#fff078", "#fff181", "#fff28b", "#fff395", "#fff6b2",
                  "#ffefaf", "#feebd2"]
    stage_chrono_names = ["Ber.", "Val.", "Hau.", "Bar.", "Apt.", "Alb.", "Cen.", "Tur.", "Con.", "San.", "Cam.",
                          "Maa.", "Dan.", "Sel.", "Tha.", "Ypr.", "Lut.", "Bart.", "Pri.", "Rup.", "Cha.", "Aqu.",
                          "Bur.", "Lan.", "Ser.", "Tor.", "Mes.", "", "", ""]

    ax["D"].bar(stage_chrono, bar_vals_stage, width=width_stage, color=cols_stage, edgecolor="none")
    ax["D"].bar_label(ax["D"].containers[0], labels=stage_chrono_names, padding=-15, fontsize = 6, rotation = 90)
    ax["D"].set_xlim(xmin= -145)
    ax["D"].set_ylim(ymin = 0, ymax = 0.4)
    ax["D"].set_yticks(ticks = [0], labels = [""])
    ax["D"].axis("off")
    # ax["D"].set_xticklabels([])

    #Time line - Epochs

    epoch_chrono = [-122.5, -83, -61, -44.95, -28.465, -14.1815, -3.9565, -1.29585, -0.00585]
    bar_vals = [0.25]*9
    width = [44.5, 34, 10, 22.1, 10.9, 17.7, 2.753, 2.5683, 0.0117]
    cols = ["#8fd07b", "#67a155", "#ffb17b", "#ffbc87", "#ffc694", "#ffeb3e", "#fff6b2", "#ffefaf", "#feebd2"]
    epoch_chrono_names = ["Early Cretaceous", "Late Cretaceous", "Paleoc.", "Eocene", "Oligocene", "Miocene", "Pli.",
                          "Ple.", ""]

    ax["E"].bar(epoch_chrono, bar_vals, width=width, color=cols, edgecolor="none")
    ax["E"].bar_label(ax["E"].containers[0], labels=epoch_chrono_names, padding=-10, fontsize = 6)
    ax["E"].set_xlim(xmin=-145)
    ax["E"].set_ylim(ymin=0, ymax=0.25)
    ax["E"].set_yticks(ticks=[0], labels=[""])

    ax["E"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
               labels=["", "", "", "", "", ""])
    ax["E"].set_xlabel("Time (Ma)", fontsize = 7)
    ax["E"].axis("off")
    ax["E"].set_visible(True)

    ax["F"].plot(time_es, freq_es, color="#bf3e36")
    ax["F"].hlines(bf2, xmin= -145, xmax = 0, linestyle = "dashed", color = "black", linewidth = 1)
    ax["F"].hlines(bf6, xmin= -150, xmax = 0, linestyle = "dashed", color = "grey", linewidth = 1)
    ax["F"].set_xlabel("Time (Ma)", fontsize=7)
    ax["F"].set_ylabel("Frequency of rate shift", fontsize=7)
    ax["F"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
               labels=["", "", "", "", "", ""])

    ax["G"].plot(time_ss, freq_ss, color="#2d7096")
    ax["G"].hlines(bf2, xmin= -145, xmax = 0, linestyle = "dashed", color = "black", linewidth = 1)
    ax["G"].hlines(bf6, xmin= -150, xmax = 0, linestyle = "dashed", color = "grey", linewidth = 1)
    ax["G"].set_xlabel("Time (Ma)", fontsize=7)
    # ax["G"].set_ylabel("Speciation\nrate shift", fontsize=7)
    ax["G"].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
               labels=["", "", "", "", "", ""])

    df = DataFrame(list(zip(time_s, rate_s, minHPD_s, maxHPD_s,  time_e, rate_e, minHPD_e, maxHPD_e)),
                      columns=['Time_s', 'Rate_s', "HPD_Min_s", "HPD_Max_s", "Time_e", "Rate_e", "HPD_Min_e", "HPD_Max_e"])


    df.to_excel("/path_to_pyrate_output_folder/{i}/rates.xlsx".format(i=file_name))

    # fig.savefig("/Users/kristinakocakova/PycharmProjects/Diversification_and_Age-dependent_Extinction/figures/all_species_{i}/rates_{i}.pdf".format(i=file_name), bbox_inches='tight')

    plt.tight_layout()
    plt.show()



# Individual rate shifts

r = read_csv("/path_to_pyrate_output_folder/{i}/RTT_plots.r".format(i = file_name), sep="\t", header= None)

# CONVERT PYRATE GENERATED R FILE INTO PYTHON LISTS OF SPECIATION AND EXTINCTION RATE VALUES
def convert(object):
    lyst = object.to_list()
    for i in lyst:
        lyst = i.split(",")
    lyst[-1] = lyst[-1].replace(")", "")
    lyst[0] = re.sub(".*" + "\(", "", lyst[0])
    for i in range(len(lyst)):
        lyst[i] = float(lyst[i])
    return lyst
def convert_bf(object):
    lyst = object.to_list()
    for i in lyst:
        x = i.split()
    y = float(x[2])
    return y

"""
IMPORTANT
The row numbers may differ depending on the output file from PyRate, which in turn depends on the parameters.
So it's best to first check in the .r file where exactly rates, CIs and times are stored. 
"""

#Extinction rate shifts
time_es = convert(r.iloc[28])
freq_es = convert(r.iloc[29])
bf2 = convert_bf(r.iloc[15])
bf6 = convert_bf(r.iloc[16])

#Speciation rate shifts
time_ss = convert(r.iloc[12])
freq_ss = convert(r.iloc[13])

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

sns.set_style("whitegrid")

fig, ax = plt.subplots(2, 1, figsize=(6.69291, 4), sharey= True)

ax[0].plot(time_es, freq_es, color="#bf3e36")
ax[0].hlines(bf2, xmin=-150, xmax=0, linestyle="dashed", color="grey", linewidth=1)
ax[0].hlines(bf6, xmin=-150, xmax=0, linestyle="dashed", color="black", linewidth=1)
ax[0].set_xlabel("Time (Ma)", fontsize=7)
ax[0].set_ylabel("Frequency of rate shift", fontsize=7)
ax[0].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])

ax[1].plot(time_ss, freq_ss, color="#2d7096")
ax[1].hlines(bf2, xmin=-150, xmax=0, linestyle="dashed", color="grey", linewidth=1)
ax[1].hlines(bf6, xmin=-150, xmax=0, linestyle="dashed", color="black", linewidth=1)
ax[1].set_xlabel("Time (Ma)", fontsize=7)
# ax["G"].set_ylabel("Speciation\nrate shift", fontsize=7)
ax[1].set_xticks(ticks=[-100.5, -66, -56, -33.9, -23.03, -5.333],
                   labels=["", "", "", "", "", ""])

legend_elements = [plt.Line2D([0], [0], color='black', lw=1, label='BF6', linestyle="dashed"), plt.Line2D([0], [1], color='grey', lw=1, label='BF2', linestyle="dashed")]
ax[0].legend(handles=legend_elements, loc="upper right", fontsize=7)

sns.despine()
plt.tight_layout()