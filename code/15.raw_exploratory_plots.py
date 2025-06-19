"""
Project: Diversification Rates and ADE
Description:
Plotting of longevity distributions based on PyRate speciation and extinction times, and on first and last appearances
Only extict taxa plotted
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# load input file

species = pd.read_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/inputs/2025/June/species.txt", sep = "\t")

# dictionary to find the first and last appearances

dyct = {}

for i in species["taxon_name"]:
    if i not in dyct:
        dyct[i] = [0, 145, 0]


for j in np.unique(species["taxon_name"]):
    subset = species.loc[species["taxon_name"] == j]
    dyct[j][2] = subset["status"].values[0]
    for k in range(len(subset)):
        if subset["max_ma"].values[k] > dyct[j][0]:
            dyct[j][0] = subset["max_ma"].values[k]

for j in np.unique(species["taxon_name"]):
    subset = species.loc[species["taxon_name"] == j]
    s = subset["status"].values[0]
    for k in range(len(subset)):
        if s == "extant":
            dyct[j][1] = 0
        elif s == "extinct":
            if subset["min_ma"].values[k] < dyct[j][1]:
                dyct[j][1] = subset["min_ma"].values[k]



species_l = list(dyct.keys())
fad = []
lad = []
status = []

for i in species_l:
    fad.append(dyct[i][0])
    lad.append(dyct[i][1])
    status.append(dyct[i][2])

df = pd.DataFrame()
df["species"] = species_l
df["fad"] = fad
df["lad"] = lad
df["status"] = status

df = df.sort_values(by = "species")

# calculate duration based on fads and lads

df["age"] = df["fad"] - df["lad"]

df.to_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species_raw.txt", sep = "\t", index = False)

# repeat for genera

# load input file

genera = pd.read_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/inputs/2025/June/genera.txt", sep = "\t")

# dictionary to find the first and last appearances

dyct_g = {}

for i in genera["taxon_name"]:
    if i not in dyct_g:
        dyct_g[i] = [0, 145, 0]


for j in np.unique(genera["taxon_name"]):
    subset = genera.loc[genera["taxon_name"] == j]
    dyct_g[j][2] = subset["genus_status"].values[0]
    for k in range(len(subset)):
        if subset["max_ma"].values[k] > dyct_g[j][0]:
            dyct_g[j][0] = subset["max_ma"].values[k]

for j in np.unique(genera["taxon_name"]):
    subset = genera.loc[genera["taxon_name"] == j]
    s = subset["genus_status"].values[0]
    for k in range(len(subset)):
        if s == "extant":
            dyct_g[j][1] = 0
        elif s == "extinct":
            if subset["min_ma"].values[k] < dyct_g[j][1]:
                dyct_g[j][1] = subset["min_ma"].values[k]



genera_l = list(dyct_g.keys())
fad = []
lad = []
status = []

for i in genera_l:
    fad.append(dyct_g[i][0])
    lad.append(dyct_g[i][1])
    status.append(dyct_g[i][2])

df_g = pd.DataFrame()
df_g["genus"] = genera_l
df_g["fad"] = fad
df_g["lad"] = lad
df_g["status"] = status

df_g = df_g.sort_values(by = "genus")

# calculate duration based on fads and lads

df_g["age"] = df_g["fad"] - df_g["lad"]

df_g.to_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/genera_raw.txt", sep = "\t", index = False)


# plot as with the PyRate ages
# only include extinct taxa
df = df.loc[df["lad"] != 0]
df_g = df_g.loc[df_g["lad"] != 0]

age_s = df["age"].values
age_g = df_g["age"].values

mean_age_s = np.mean(age_s)
mean_age_g = np.mean(age_g)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6
fig, ax = plt.subplots(1 , 2, figsize = (6.69291, 8.26772/3))
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)


sns.set_style("whitegrid")
h = ax[0].hist(age_s, edgecolor = "white", color = "#a8293c", bins = 20)
ymax = np.max(h[0])
sns.despine()
ax[0].vlines(x=[mean_age_s], ymin=0, ymax = ymax, linestyles="dashed", colors=['#000000'], alpha=0.5, linewidth = 2)
ax[0].text(mean_age_s + 3, ymax - 15, s =str(round(mean_age_s, 3)), c = '#000000', fontsize = 9)

h_g = ax[1].hist(age_g, edgecolor = "white", color = "#a8293c", bins = 20)
ymax_g = np.max(h_g[1])
sns.despine()
ax[1].vlines(x=[mean_age_g], ymin=0, ymax = ymax_g, linestyles="dashed", colors=['#000000'], alpha=0.5, linewidth = 2)
ax[1].text(mean_age_g + 3, ymax_g - 10, s =str(round(mean_age_g, 3)), c = '#000000', fontsize = 9)


fig.supxlabel("Longevity (Myr)", fontsize=10)
fig.supylabel("Frequency", fontsize=10)

plt.tight_layout()
plt.show(block = True)





