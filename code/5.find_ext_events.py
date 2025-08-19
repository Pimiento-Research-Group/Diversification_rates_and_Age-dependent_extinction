"""
Project: Diversification Rates and ADE
Description:
Calculate the background rate and identify times when the mean rate reaches at least 3 or 6 times the background
"""

from pandas import *
import re
import numpy as np
np.set_printoptions(suppress=True, precision=3)

f = "path_to_rates.xlsx" # "/Volumes/External_memory/Dropbox/Kristina_PhD_K_version/Kristina_files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species"

rates = read_excel(f + "/rates.xlsx")

rates["Time_e"] = rates["Time_e"] * -1
rates["Time_s"] = rates["Time_s"] * -1
rates["Time_d"] = rates["Time_d"] * -1

time_ext = np.array(rates["Time_e"])
rate_ext = np.array(rates["Rate_e"])
time_sp = np.array(rates["Time_s"])
rate_sp = np.array(rates["Rate_s"])
time_d = np.array(rates["Time_d"])
rate_d = np.array(rates["Rate_d"])

harm_mean_e = len(rate_ext)/np.sum(1/np.array(rate_ext))
harm_mean_s = len(rate_sp)/np.sum(1/np.array(rate_sp))

time_ext[rate_ext > 3*harm_mean_e] # 6* for intense ext events
time_sp[rate_sp > 3*harm_mean_s] # 6* for intense ext events
time_d[rate_d > 3* 0.02]
time_d[rate_d < 3* -0.02]

ext = [[145, 95.71], [95.71, 83.31], [83.31, 73.40], [73.40, 71.90], [71.90, 66.60], [66.60, 65.79],
              [65.79, 55.89], [55.89, 38.58], [38.58, 33.57], [33.57, 3.65], [3.65, 0.0117]]

spec = [[145, 112.23], [112.23, 101.22], [101.22, 100.12], [100.12, 87.21], [87.21, 82.70], [82.70, 72.40],
        [72.40, 71.39], [71.39, 66.17], [66.17, 64.69], [64.69, 56.49], [56.49, 55.59]]

div = [[145, 111.34], [111.34, 101.42], [101.42, 99.72], [99.72, 94.42], [94.42, 93.62], [93.62, 87.11], [87.11, 82.61],
       [82.61, 73.60], [73.60, 72.20], [72.20, 72.10], [72.10, 70.80], [70.80, 66.69], [66.69, 65.79], [65.79, 65.59],
       [65.59, 64.99], [64.99, 56.49], [56.49, 55.49], [55.49, 38.28], [38.28, 33.37], [33.37, 3.15], [3.15, 0.0117]]

# maximum rates + CIs (Table 1)

for i in ext:
    sub_rate = rate_ext[(time_ext <= i[0]) & (time_ext > i[1])]
    m = np.max(sub_rate)
    r = rates[["Time_e", "Rate_e", "HPD_Min_e", "HPD_Max_e"]]
    print(i)
    print(r[r["Rate_e"] == m])

for i in spec:
    sub_rate = rate_sp[(time_sp <= i[0]) & (time_sp > i[1])]
    m = np.max(sub_rate)
    r = rates[["Time_s", "Rate_s", "HPD_Min_s", "HPD_Max_s"]]
    print(i)
    print(r[r["Rate_s"] == m])

for i in div:
    sub_rate = rate_d[(time_d <= i[0]) & (time_d > i[1])]
    max = np.max(sub_rate)
    min = np.min(sub_rate)
    r = rates[["Time_d", "Rate_d", "HPD_Min_d", "HPD_Max_d"]]
    print(i)
    print(r[r["Rate_d"] == max])
    print(r[r["Rate_d"] == min])




