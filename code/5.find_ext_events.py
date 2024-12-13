"""
Project: Diversification Rates and ADE
Author: Kristína Kocáková
Description:
Calculate the background rate and identify times when the mean rate reaches at least 3 or 6 times the background
"""

from pandas import *
import re
import numpy as np
np.set_printoptions(suppress=True, precision=3)

def extract_ext_rate(directory):
    r = read_csv(
        "/Users/kristinakocakova/PycharmProjects/Diversification_and_Age-dependent_Extinction/data/all_species{i}/RTT_plots.r".format(
            i=directory), sep="\t", header=None)

    def convert(object):
        lyst = object.to_list()
        for i in lyst:
            lyst = i.split(",")
        lyst[-1] = lyst[-1].replace(")", "")
        lyst[0] = re.sub(".*" + "\(", "", lyst[0])
        for i in range(len(lyst)):
            lyst[i] = float(lyst[i])
        return lyst

    time_e = convert(r.iloc[19]) #for origination channge to 3
    rate_e = convert(r.iloc[20]) #for origination change to 4

    return time_e, rate_e

time_ext, rate_ext = extract_ext_rate("")

time_ext = np.array(time_ext)
rate_ext = np.array(rate_ext)

harm_mean = len(rate_ext)/np.sum(1/np.array(rate_ext))

time_ext[rate_ext > 3*harm_mean]









