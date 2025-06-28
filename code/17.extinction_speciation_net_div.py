"""
Project: Diversification Rates and ADE
Description:
Calculate extinction, speciation and net diversification based on FADs and LADs
https://www.jstor.org/stable/1571654?seq=3
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import statsmodels.api as sm
import statsmodels.formula.api as smf
from numpy.ma.extras import row_stack
from collections import Counter
from scipy.stats import skew

species = pd.read_csv("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species_raw.txt", sep = "\t")

time_ints = [145, 139.8, 132.6, 125.77, 121.4, 113, 100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66, 61.6, 59.2, 56, 47.8, 41.2, 37.71, 33.9, 27.82, 23.03, 20.44, 15.98, 13.82, 11.63, 7.246,
                5.333, 3.6, 2.58, 1.8, 0.774, 0.129, 0.0117, 0]

for i in range(len(time_ints) - 1):
    df = species.loc[(species["lad"] <= time_ints[i]) & (species["lad"] >= time_ints[i + 1])]


