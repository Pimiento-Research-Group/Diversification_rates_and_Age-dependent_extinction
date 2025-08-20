"""
Project: Diversification Rates and ADE
Description:
Prediction of the Weibull distribution shape parameter based on empirical data using the trained models created in script 6.
"""

import glob, os, argparse, sys
import numpy as np
import pandas as pd
import scipy.special
import scipy.ndimage as nd
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from scipy import stats
from sklearn.metrics import confusion_matrix

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import metrics
import matplotlib.pyplot as plt
import pickle as pkl
import statistics
import itertools
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
np.set_printoptions(suppress=True, precision=3)  # prints floats, no scientific notation


def print_update(s):
    sys.stdout.write('\r')
    sys.stdout.write(s)
    sys.stdout.flush()


def print_update(s):
    sys.stdout.write('\r')
    sys.stdout.write(s)
    sys.stdout.flush()


def data_simulator_mixture(q=None,  # if None: rnd drawn
                           N=1000,
                           min_data_size=1,
                           max_data_size=1000,
                           fixed_shape=None,  # if not None -> array of 2 values
                           fixed_scale=None,  # if not None -> array of 2 values
                           magnitude=2,  # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                           min_longevity=2,  # ~ exp( U[log(min_longevity), log(max_longevity)] )
                           max_longevity=30,  # ~ exp( U[log(min_longevity), log(max_longevity)] )
                           n_Hbins=100,  # bins of histogram occs-per-species
                           maxNbinsPerSpecies=20,  # bins of histogram timebins-per-species
                           gamma_model=True,
                           alpha=1,  # rate heterogeneity across species (smaller alpha = greater variation)
                           verbose=0,
                           plot=False,
                           include_info=True,
                           single_longevity=True,
                           single_shape=True,
                           two_shapes=False,
                           mixture=False):
    """
    N: number of datasets


    """
    Hbins = np.linspace(0.5, (n_Hbins + 0.5), n_Hbins + 1)
    Hbins = np.array(list(Hbins) + [1000.5])
    EpochBinSize = np.random.uniform(2, 10, N)

    # MODEL SETTINGS
    print("\n")

    # preservation
    if q is None:
        q_rates = np.random.uniform(0.25, 1, N)
    else:
        q_rates = np.ones(N) * q

    hSize = len(Hbins) - 1 + maxNbinsPerSpecies + 2
    features = np.zeros((N, hSize))
    n_species = np.random.randint(min_data_size, max_data_size, N)

    if fixed_shape is None:
        rnd_shapes = np.zeros((N, 2))
        if single_shape == True:
            for i in range(N):
                rnd_shape = np.exp(np.random.uniform(np.log(1 / magnitude), np.log(magnitude)))
                rnd_shape = np.tile(rnd_shape, (1, 2))
                rnd_shape = np.sort(rnd_shape, 1)
                rnd_shapes[i] = rnd_shape
        elif two_shapes == True:
            for i in range(N):
                rnd_shape_1 = np.exp(np.random.uniform(np.log(1 / magnitude), 0))
                rnd_shape_2 = np.exp(np.random.uniform(0, np.log(magnitude)))
                rnd_shapes[i] = np.array([rnd_shape_1, rnd_shape_2])
        elif mixture == True:
            for i in range(N):
                x = np.random.randint(0, 3)
                if x == 0:  # both shapes == 1
                    rnd_shape = np.ones((1))
                    rnd_shape = np.tile(rnd_shape, (1, 2))
                    rnd_shapes[i] = rnd_shape
                elif x == 1:  # both shapes identical
                    rnd_shape = np.exp(np.random.uniform(np.log(1 / magnitude), np.log(magnitude)))
                    rnd_shape = np.tile(rnd_shape, (1, 2))
                    rnd_shape = np.sort(rnd_shape, 1)
                    rnd_shapes[i] = rnd_shape
                elif x == 2:  # two different shapes
                    rnd_shape_1 = np.exp(np.random.uniform(np.log(1 / magnitude), 0))
                    rnd_shape_2 = np.exp(np.random.uniform(0, np.log(magnitude)))
                    rnd_shapes[i] = np.array([rnd_shape_1, rnd_shape_2])

    else:
        if single_shape:
            rnd_shapes = np.ones((N, 1)) * fixed_shape
            rnd_shapes = np.tile(rnd_shapes, (1, 2))
        else:
            rnd_shapes = np.ones((N, 2)) * fixed_shape

    if fixed_scale is None:
        if single_longevity:
            mean_longevity = np.exp(np.random.uniform(np.log(min_longevity), np.log(max_longevity), (N, 1)))
        else:
            # draw rnd scale
            mean_longevity = np.exp(np.random.uniform(np.log(min_longevity), np.log(max_longevity), (N, 2)))

        # mean_longevity = np.sort(mean_longevity, 1)
        rnd_scales = mean_longevity / scipy.special.gamma(1 + 1 / rnd_shapes)
    else:
        rnd_scales = np.ones((N, 2)) * fixed_scale
        mean_longevity = rnd_scales * scipy.special.gamma(1 + 1 / rnd_shapes)

    nFsampled = []
    for i in np.arange(N):
        # loop over datasets
        EpochBins = np.linspace(0, 1000, 1 + int(1000 / EpochBinSize[i]))
        # if i % 10 == 0:
        #     print_update("Simulating data...%s/%s" % (i, N))

        # simulate species longevities
        compound_index = np.random.binomial(1, 0.5, n_species[i])  # assign species to one of the two Weibulls
        rW = np.random.weibull(rnd_shapes[i][compound_index], n_species[i]) * rnd_scales[i][compound_index]
        rW[rW > 990] = 990  # avoid extreme lifespans falling outside of the allowed time window
        if plot:
            plt.hist(rW)
            plt.show()

        # simulate fossil record (n. occs)
        q_dataset = q_rates[i]
        if gamma_model:
            q_dataset = np.random.gamma(alpha, 1 / alpha, len(rW))

        # actual occs| hist of how many sp extend over 1, 2, ... time bins
        brL = 0
        while brL == 0:
            nF = np.random.poisson(q_rates[i] * rW)
            rnd_spec_time = np.random.uniform(0, 10, n_species[i])
            hOccs = np.zeros(maxNbinsPerSpecies)
            for j in range(n_species[i]):
                if nF[j] > 0:
                    # species with no records are not obs data
                    oF = rnd_spec_time[j] + np.random.uniform(0, rW[j], nF[j])
                    ndigi = np.digitize(oF, EpochBins)
                    nBins = np.min([maxNbinsPerSpecies, len(np.unique(ndigi))]) - 1
                    hOccs[nBins] = hOccs[nBins] + 1
                    coarse_occs = np.random.uniform(EpochBins[ndigi - 1], EpochBins[ndigi])
                    br_temp = np.max(coarse_occs) - np.min(coarse_occs)
                    brL += br_temp

            if brL == 0:
                print(q_rates[i])
                q_estimate = None
                # sys.exit("???")
            else:
                q_estimate = np.sum(nF[nF > 1] - 1) / brL
        if verbose:
            print(q_rates[i], q_estimate)

        # BUILDING FEATURE SET
        # fraction of species with occurrences in 1, 2, 3, ... N time bins
        hOccs = hOccs / np.sum(hOccs)

        nF[nF > 1000] = 1000
        nFsampled.append(len(nF[nF > 0]))
        if maxNbinsPerSpecies > 0:
            # append fraction of species with 1, 2, 3 ... N occurrences
            h_temp = np.append(np.histogram(nF[nF > 0], bins=Hbins, density=True)[0], hOccs)
            # append avg bin size
            h_temp = np.append(h_temp, EpochBinSize[i])
            # append approximate preservation rate
            features[i, :] = np.append(h_temp, q_estimate)
        else:
            h_temp = np.histogram(nFsampled, bins=Hbins, density=True)[0]
            features[i, :] = np.append(h_temp, EpochBinSize[i])

    # BUILD LABELS
    labels = np.hstack((rnd_shapes, mean_longevity))
    print("\nDone.")
    if include_info:
        info = {"EpochBinSize": EpochBinSize,
                "n_species": n_species,
                "q_rates": q_rates,
                "rnd_scales": rnd_scales,
                "mean_longevity": mean_longevity,
                "sample_species": np.array(nFsampled)
                }
    else:
        info = None
    return features, labels, info


def calcHPD(data, level):
    assert (0 < level < 1)
    d = list(data)
    d.sort()
    nData = len(data)
    nIn = int(round(level * nData))
    if nIn < 2:
        sys.exit('\n\nToo little data to calculate marginal parameters.')
    i = 0
    r = d[i + nIn - 1] - d[i]
    for k in range(len(d) - (nIn - 1)):
        rk = d[k + nIn - 1] - d[k]
        if rk < r:
            r = rk
            i = k
    assert 0 <= i <= i + nIn - 1 < len(d)
    return (d[i], d[i + nIn - 1])


def get_mean_bin_size(tbl):
    tbl_sp = tbl[1:, 2:4].astype(float)
    oF = np.max(tbl_sp, 1) - np.min(tbl_sp, 1)
    return np.mean(oF)


def estimate_q_from_range_data(tbl):
    minAges = np.min(tbl[:, 2:4].astype(float), axis=1)
    maxAges = np.max(tbl[:, 2:4].astype(float), axis=1)
    # assign species indexes
    sp_id = tbl[:, 0]
    sp_id_num = pd.factorize(sp_id)[0]

    # count occs per species
    num_occs_per_species = np.unique(sp_id_num, return_counts=1)[1]

    # get branch length
    ages = np.random.uniform(minAges, maxAges)
    m1 = nd.maximum(ages, sp_id_num, np.unique(sp_id_num))
    m2 = nd.minimum(ages, sp_id_num, np.unique(sp_id_num))
    rndBr = m1 - m2
    qA = np.sum(num_occs_per_species[num_occs_per_species > 1] - 1) / np.sum(rndBr[num_occs_per_species > 1] + 0.00001)
    # MLE of q rate
    qA = np.sum(num_occs_per_species[num_occs_per_species > 1] - 1) / np.sum(rndBr[num_occs_per_species > 1] + 0.00001)
    return qA


def get_data_array(f, time_slice=None, n_rnd_q_estimates=1, n_Hbins=100, maxNbinsPerSpecies=20, min_n_taxa=1):
    Hbins = np.linspace(0.5, (n_Hbins + 0.5), n_Hbins + 1)
    Hbins = np.array(list(Hbins) + [1000.5])

    f_name = os.path.splitext(os.path.basename(f))[0]
    rnd_data = []
    for i in range(n_rnd_q_estimates):
        tbl = pd.read_csv(f, sep="\t").to_numpy()

        if time_slice is not None:
            n_remaining_sp, remaining_sp_names, tbl = filter_species_time_slice(tbl, time_slice)
            if n_remaining_sp < min_n_taxa:
                print("Not enough species in time slice (%s)" % n_remaining_sp)
                return [0, 0, 0, 0, 0]
            elif i == 0:
                print(n_remaining_sp, "species remaining after filtering time slice")

        sp_label = np.unique(tbl[1:, 0])
        num_occs_per_species = np.unique(tbl[1:, 0], return_counts=1)[1]
        hist = np.histogram(num_occs_per_species, bins=Hbins, density=True)[0]

        bin_size = get_mean_bin_size(tbl)
        EpochBins = np.linspace(0, 1000, 1 + int(1000 / bin_size))

        # INTRODUCE TAXONOMIC BIAS
        if False:
            num_occs_per_species_pr = 1 / num_occs_per_species
            pr_synonym = num_occs_per_species_pr / np.sum(num_occs_per_species_pr)

            # synomize
            n_synonyms = int(run_empirical_taxonbias * len(sp_label))
            syn_sp_lab = np.random.choice(sp_label, n_synonyms, p=pr_synonym,
                                          replace=False)  # species to be turned into new label
            other_sp_lab = np.setdiff1d(sp_label, syn_sp_lab)
            new_names = np.random.choice(other_sp_lab, n_synonyms, replace=True)
            tax_i = 0
            for tax in syn_sp_lab:
                tbl[tbl[:, 0] == tax, 0] = new_names[tax_i]

        sp_label = np.unique(tbl[1:, 0])
        num_occs_per_species = np.unique(tbl[1:, 0], return_counts=1)[1]
        hist = np.histogram(num_occs_per_species, bins=Hbins, density=True)[0]

        q_est = estimate_q_from_range_data(tbl)
        # count boundary crossers
        hOccs = np.zeros(maxNbinsPerSpecies)
        for tax in sp_label:
            tbl_sp = tbl[tbl[:, 0] == tax, 2:4].astype(float)
            # randomize fossil age
            oF = np.random.uniform(np.min(tbl_sp, 1), np.max(tbl_sp, 1))
            nBins = np.min([maxNbinsPerSpecies, len(np.unique(np.digitize(oF, EpochBins)))]) - 1
            hOccs[nBins] = hOccs[nBins] + 1

        hOccs = hOccs / np.sum(hOccs)
        # print(hOccs)
        hist = np.append(hist, hOccs)
        hist_temp = np.append(hist, bin_size)

        hist = np.append(hist_temp, q_est)
        rnd_data.append(hist)

    rnd_data = np.array(rnd_data)
    [binSize, Q_est] = np.mean(rnd_data, axis=0)[-2:]
    output = {'features': rnd_data,
              'n_species': len(sp_label),
              'species_names': remaining_sp_names,
              'n_occs': np.sum(num_occs_per_species),
              'bin_size': binSize,
              'estimated_q': Q_est
              }
    return output


class AdeNNrescaler():
    def __init__(self,
                 log_shapes=True,
                 log_longevities=True,
                 shapes_mul=1,
                 longevities_mul=1,
                 shape_ratio=False
                 ):

        self.log_shapes = log_shapes
        self.log_longevities = log_longevities
        self.shapes_mul = shapes_mul
        self.longevities_mul = longevities_mul
        self.shape_ratio = shape_ratio

    def transform_labels(self, shapes, longevities, stack=False):
        tr_shapes = shapes + 0
        tr_longevities = longevities + 0

        if self.shape_ratio:
            tr_shapes[:, 1] = shapes[:, 1] / shapes[:, 0]

            if self.log_shapes:
                tr_shapes[:, 0] = np.log(tr_shapes[:, 0])

        elif self.log_shapes:
            tr_shapes = np.log(tr_shapes)

        if self.log_longevities:
            tr_longevities = np.log(tr_longevities)

        tr_longevities *= self.longevities_mul
        tr_shapes *= self.shapes_mul

        if stack:
            lab = np.hstack((tr_shapes, tr_longevities))
            return lab
        else:
            return tr_shapes, tr_longevities

    def back_transform_labels(self, tr_shapes, tr_longevities, stack=False):
        longevities = tr_longevities / self.longevities_mul
        shapes = tr_shapes / self.shapes_mul

        if self.shape_ratio:
            if self.log_shapes:
                shapes[:, 0] = np.exp(shapes[:, 0])

            shapes[:, 1] = tr_shapes[:, 1] * shapes[:, 0]

        elif self.log_shapes:
            shapes = np.exp(shapes)

        if self.log_longevities:
            longevities = np.exp(longevities)

        if stack:
            lab = np.hstack((shapes, longevities))
            return lab
        else:
            return shapes, longevities


def filter_species_time_slice(tbl, time_slice):
    max_occs_ages = np.max(tbl[:, 2:4].astype(float), 1)
    min_occs_ages = np.min(tbl[:, 2:4].astype(float), 1)
    sp_name, sp_indx = np.unique(tbl[:, 0], return_inverse=True)
    sp_include = []
    n_species_included = []
    for i in np.unique(sp_indx):
        max_occs_ages_sp = max_occs_ages[sp_indx == i]
        min_occs_ages_sp = min_occs_ages[sp_indx == i]
        sp_ext_time = np.min(min_occs_ages_sp)
        if sp_ext_time <= np.max(time_slice) and sp_ext_time >= np.min(time_slice):
            sp_include = sp_include + list(np.where(sp_indx == i)[0])
            n_species_included.append(sp_name[i])
    sp_include_indx_in_tbl = np.array(sp_include)

    tbl_new = tbl[sp_include_indx_in_tbl]
    return len(np.unique(n_species_included)), sp_include, tbl_new


def labels_to_class(labels_original, rescaler=None):
    class_labels = np.zeros((labels_original.shape[0], 3))
    labels = labels_original + 0
    if rescaler is not None:
        if rescaler.log_shapes:
            labels[:, 0] = np.exp(labels[:, 0])
            labels[:, 1] = np.exp(labels[:, 1])

    class_labels[np.where((labels[:, 0] == 1) & (labels[:, 1] == 1)), 0] = 1  # no ADE
    class_labels[np.where((labels[:, 0] == labels[:, 1]) & (labels[:, 0] != 1)), 1] = 1  # single ADE
    class_labels[np.where((np.sum(class_labels, 1) == 0)), 2] = 1  # double ADE

    return class_labels


def calcCI(predicted_labels, errors, mse=False):
    if mse == True:
        errors = np.sqrt(errors)
        CI1 = (predicted_labels - 1.96 * errors)
        CI2 = (predicted_labels + 1.96 * errors)
    else:  # i.e. if errors are rmse
        CI1 = (predicted_labels - 1.96 * errors)
        CI2 = (predicted_labels + 1.96 * errors)
    return np.array([CI1, CI2])


def get_rmse(taxa_no, model_type):  # model type - 1, 2
    if model_type == 1:
        n_pred = 1000
        features_arr = np.zeros((1, 123))
        labels_arr = np.zeros((1, 3))
        taxa_10p = (taxa_no/10)
        while features_arr.shape[0] < n_pred:
            features, labels, info = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                            N=5000,
                                                            min_data_size=(taxa_no + taxa_10p),
                                                            max_data_size=(taxa_no + taxa_10p*2),
                                                            fixed_shape=None,  # if not None -> array of 2 values
                                                            fixed_scale=None,  # if not None -> array of 2 values
                                                            magnitude=2,
                                                            # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                            min_longevity=2,
                                                            # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                            max_longevity=30,
                                                            # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                            n_Hbins=100,  # bins of histogram occs-per-species
                                                            maxNbinsPerSpecies=20,
                                                            # bins of histogram timebins-per-species
                                                            gamma_model=True,
                                                            alpha=1,
                                                            # rate heterogeneity across species (smaller alpha = greater variation)
                                                            verbose=0,
                                                            plot=False,
                                                            include_info=True,
                                                            single_longevity=True,
                                                            single_shape=True,
                                                            two_shapes=False,
                                                            mixture=False)
            x = info["sample_species"]
            indx = np.where((x > taxa_no - taxa_10p) & (x < taxa_no + taxa_10p))[0]
            features = features[indx, :]
            labels = labels[indx, :]
            features_arr= np.append(features_arr, features, axis = 0)
            labels_arr= np.append(labels_arr, labels, axis = 0)
            print(features_arr.shape[0])

            if features_arr.shape[0] >= n_pred:
                features_arr = features_arr[1:, :]
                labels_arr = labels_arr[1:, :]
                print("ADE1 " + str(features_arr.shape[0]))
                break


        rescaler = AdeNNrescaler(log_shapes=True,
                                 log_longevities=True,
                                 shapes_mul=1,
                                 longevities_mul=1,
                                 shape_ratio=False)
        shapes_r, long_r = rescaler.transform_labels(labels_arr[:, :2], labels_arr[:, 2:])
        labels_r = np.hstack((shapes_r, long_r))

        n_predictions = 1000  #
        pred = np.array([model_single_shape(features_arr, training=True) for _ in range(n_predictions)])
        pred_mean = np.mean(pred, axis=0)

        rmse = np.sqrt(np.mean((labels_r - pred_mean) ** 2, axis=0))

    elif model_type == 2:
        n_pred = 1000
        features_arr = np.zeros((1, 123))
        labels_arr = np.zeros((1, 3))
        taxa_10p = (taxa_no/10)
        while features_arr.shape[0] < n_pred:
            features, labels, info = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                            N=5000,
                                                            min_data_size=(taxa_no + taxa_10p),
                                                            max_data_size=(taxa_no + taxa_10p*2),
                                                            fixed_shape=None,  # if not None -> array of 2 values
                                                            fixed_scale=None,  # if not None -> array of 2 values
                                                            magnitude=2,
                                                            # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                            min_longevity=2,
                                                            # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                            max_longevity=30,
                                                            # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                            n_Hbins=100,  # bins of histogram occs-per-species
                                                            maxNbinsPerSpecies=20,
                                                            # bins of histogram timebins-per-species
                                                            gamma_model=True,
                                                            alpha=1,
                                                            # rate heterogeneity across species (smaller alpha = greater variation)
                                                            verbose=0,
                                                            plot=False,
                                                            include_info=True,
                                                            single_longevity=True,
                                                            single_shape=False,
                                                            two_shapes=True,
                                                            mixture=False)
            x = info["sample_species"]
            indx = np.where((x > taxa_no - taxa_10p) & (x < taxa_no + taxa_10p))[0]
            features = features[indx, :]
            labels = labels[indx, :]
            features_arr= np.append(features_arr, features, axis = 0)
            labels_arr= np.append(labels_arr, labels, axis = 0)
            print(features_arr.shape[0])

            if features_arr.shape[0] >= n_pred:
                features_arr = features_arr[1:, :]
                labels_arr = labels_arr[1:, :]
                print("ADE2 " + str(features_arr.shape[0]))
                break


        rescaler = AdeNNrescaler(log_shapes=True,
                                 log_longevities=True,
                                 shapes_mul=1,
                                 longevities_mul=1,
                                 shape_ratio=False)
        shapes_r, long_r = rescaler.transform_labels(labels_arr[:, :2], labels_arr[:, 2:])
        labels_r = np.hstack((shapes_r, long_r))

        n_predictions = 1000  #
        pred = np.array([model_two_shapes(features_arr, training=True) for _ in range(n_predictions)])
        pred_mean = np.mean(pred, axis=0)

        rmse = np.sqrt(np.mean((labels_r - pred_mean) ** 2, axis=0))

    elif model_type == 0:
        n_pred = 1000
        features_arr = np.zeros((1, 123))
        labels_arr = np.zeros((1, 3))
        taxa_10p = (taxa_no/10)
        while features_arr.shape[0] < n_pred:
            features, labels, info = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                            N=5000,
                                                            min_data_size=(taxa_no + taxa_10p),
                                                            max_data_size=(taxa_no + taxa_10p*2),
                                                            fixed_shape=1,  # if not None -> array of 2 values
                                                            fixed_scale=None,  # if not None -> array of 2 values
                                                            magnitude=2,
                                                            # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                            min_longevity=2,
                                                            # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                            max_longevity=30,
                                                            # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                            n_Hbins=100,  # bins of histogram occs-per-species
                                                            maxNbinsPerSpecies=20,
                                                            # bins of histogram timebins-per-species
                                                            gamma_model=True,
                                                            alpha=1,
                                                            # rate heterogeneity across species (smaller alpha = greater variation)
                                                            verbose=0,
                                                            plot=False,
                                                            include_info=True,
                                                            single_longevity=True,
                                                            single_shape=True,
                                                            two_shapes=False,
                                                            mixture=False)
            x = info["sample_species"]
            indx = np.where((x > taxa_no - taxa_10p) & (x < taxa_no + taxa_10p))[0]
            features = features[indx, :]
            labels = labels[indx, :]
            features_arr= np.append(features_arr, features, axis = 0)
            labels_arr= np.append(labels_arr, labels, axis = 0)
            print(features_arr.shape[0])

            if features_arr.shape[0] >= n_pred:
                features_arr = features_arr[1:, :]
                labels_arr = labels_arr[1:, :]
                print("ADE0 " + str(features_arr.shape[0]))
                break

        rescaler = AdeNNrescaler(log_shapes=True,
                                 log_longevities=True,
                                 shapes_mul=1,
                                 longevities_mul=1,
                                 shape_ratio=False)
        shapes_r, long_r = rescaler.transform_labels(labels_arr[:, :2], labels_arr[:, 2:])
        labels_r = np.hstack((shapes_r, long_r))

        n_predictions = 1000  #
        pred = np.array([model_shape_1(features_arr, training=True) for _ in range(n_predictions)])
        pred_mean = np.mean(pred, axis=0)

        rmse = np.sqrt(np.mean((labels_r - pred_mean) ** 2, axis=0))

    return rmse


def pipeline(f, n_predictions, time_slice=None):
    n_predictions = n_predictions
    features = np.zeros((n_predictions, 123))
    taxa_no = 0
    for i in range(n_predictions):
        features_dict = get_data_array(f, time_slice)
        features[i, :] = features_dict["features"]
        taxa_no = features_dict["n_species"]

    clas = model_clas.predict(features)
    clas = np.argmax(clas, axis=1)
    # select the class with the highest frequency
    count = Counter(clas)
    clas_max = max(count, key=count.get)
    clas = clas_max

    # predict shapes and longevities
    shapes_0 = np.array([model_shape_1(features, training=True)])
    shapes_0 = shapes_0.squeeze(axis=0)
    shapes_1 = np.array([model_single_shape(features, training=True)])
    shapes_1 = shapes_1.squeeze(axis=0)
    shapes_2 = np.array([model_two_shapes(features, training=True)])
    shapes_2 = shapes_2.squeeze(axis=0)

    # obtain rmse
    rmse_ade1_ = get_rmse(taxa_no, 1)
    rmse_ade1 = rmse_ade1_[0]
    rmse_ade1_l = rmse_ade1_[2]
    rmse_ade2 = get_rmse(taxa_no, 2)  # returns array with 3 values
    rmse_ade2_1 = rmse_ade2[0]  # select the first rmse_ade2 value
    rmse_ade2_2 = rmse_ade2[1]
    rmse_ade2_l = rmse_ade2[2]
    rmse_no_ade_ = get_rmse(taxa_no, 0)
    rmse_no_ade = rmse_no_ade_[0]

    rescaler = AdeNNrescaler(log_shapes=True,
                             log_longevities=True,
                             shapes_mul=1,
                             longevities_mul=1,
                             shape_ratio=False)

    # CIs for ADE1
    ci_1 = calcCI(shapes_1, errors=rmse_ade1)
    min_ci_1 = rescaler.back_transform_labels(ci_1[0][:, :2],
                                              ci_1[0][:, 2:],
                                              stack=True)
    max_ci_1 = rescaler.back_transform_labels(ci_1[1][:, :2],
                                              ci_1[1][:, 2:],
                                              stack=True)
    ci_backscaled_1 = np.array([min_ci_1, max_ci_1])

    ci_1_l = calcCI(shapes_1, errors=rmse_ade1_[2])
    min_ci_1_l = rescaler.back_transform_labels(ci_1_l[0][:, :2],
                                                ci_1_l[0][:, 2:],
                                                stack=True)
    max_ci_1_l = rescaler.back_transform_labels(ci_1_l[1][:, :2],
                                                ci_1_l[1][:, 2:],
                                                stack=True)
    ci_backscaled_1_l = np.array([min_ci_1_l, max_ci_1_l])

    # CIs for ADE2
    ci_2_1 = calcCI(shapes_2, errors=rmse_ade2_1)
    min_ci_2_1 = rescaler.back_transform_labels(ci_2_1[0][:, :2],
                                                ci_2_1[0][:, 2:],
                                                stack=True)
    max_ci_2_1 = rescaler.back_transform_labels(ci_2_1[1][:, :2],
                                                ci_2_1[1][:, 2:],
                                                stack=True)
    ci_backscaled_2_1 = np.array([min_ci_2_1, max_ci_2_1])

    ci_2_2 = calcCI(shapes_2, errors=rmse_ade2_2)
    min_ci_2_2 = rescaler.back_transform_labels(ci_2_2[0][:, :2],
                                                ci_2_2[0][:, 2:],
                                                stack=True)
    max_ci_2_2 = rescaler.back_transform_labels(ci_2_2[1][:, :2],
                                                ci_2_2[1][:, 2:],
                                                stack=True)
    ci_backscaled_2_2 = np.array([min_ci_2_2, max_ci_2_2])

    ci_2_l = calcCI(shapes_2, errors=rmse_ade2[2])
    min_ci_2_l = rescaler.back_transform_labels(ci_2_l[0][:, :2],
                                                ci_2_l[0][:, 2:],
                                                stack=True)
    max_ci_2_l = rescaler.back_transform_labels(ci_2_l[1][:, :2],
                                                ci_2_l[1][:, 2:],
                                                stack=True)
    ci_backscaled_2_l = np.array([min_ci_2_l, max_ci_2_l])

    # CIs for ADE0
    ci_0 = calcCI(shapes_0, errors=rmse_no_ade)
    min_ci_0 = rescaler.back_transform_labels(ci_0[0][:, :2],
                                              ci_0[0][:, 2:],
                                              stack=True)
    max_ci_0 = rescaler.back_transform_labels(ci_0[1][:, :2],
                                              ci_0[1][:, 2:],
                                              stack=True)
    ci_backscaled_0 = np.array([min_ci_0, max_ci_0])

    ci_0_l = calcCI(shapes_1, errors=rmse_no_ade_[2])
    min_ci_0_l = rescaler.back_transform_labels(ci_0_l[0][:, :2],
                                                ci_0_l[0][:, 2:],
                                                stack=True)
    max_ci_0_l = rescaler.back_transform_labels(ci_0_l[1][:, :2],
                                                ci_0_l[1][:, 2:],
                                                stack=True)
    ci_backscaled_0_l = np.array([min_ci_0_l, max_ci_0_l])

    # transform
    shape_r_0, long_r_0 = rescaler.back_transform_labels(shapes_0[:, :2], shapes_0[:, 2:])
    labels_shape_r_0 = np.hstack((shape_r_0, long_r_0))
    shape_r_1, long_r_1 = rescaler.back_transform_labels(shapes_1[:, :2], shapes_1[:, 2:])
    labels_shape_r_1 = np.hstack((shape_r_1, long_r_1))
    shape_r_2, long_r_2 = rescaler.back_transform_labels(shapes_2[:, :2], shapes_2[:, 2:])
    labels_shape_r_2 = np.hstack((shape_r_2, long_r_2))

    # build the final array containg all the info
    final = np.zeros((n_predictions, 21))
    # final[:, 0] = time_slice - can't add this at this point, will add it when creating a dataframe for plotting
    final[:, 1] = clas
    final[:, 2] = labels_shape_r_0[:, 2]  # longevity under no ADE
    final[:, 3] = labels_shape_r_1[:, 0]  # shapes under ADE1
    final[:, 4] = ci_backscaled_1[0, :, 0]  # min CI ADE1
    final[:, 5] = ci_backscaled_1[1, :, 0]  # max CI ADE1
    final[:, 6] = labels_shape_r_1[:, 2]  # longevity ADE1
    final[:, 7] = ci_backscaled_1_l[0, :, 2]  # min CI of longevity ADE1
    final[:, 8] = ci_backscaled_1_l[1, :, 2]  # max CI of longevity ADE1
    final[:, 9] = labels_shape_r_2[:, 0]  # shape 1 of ADE2
    final[:, 10] = ci_backscaled_2_1[0, :, 0]  # min CI of shape 1 of ADE2
    final[:, 11] = ci_backscaled_2_1[1, :, 0]  # max CI of shape 1 of ADE2
    final[:, 12] = labels_shape_r_2[:, 1]  # shape 2 of ADE2
    final[:, 13] = ci_backscaled_2_2[0, :, 1]  # min CI of shape 2 of ADE2
    final[:, 14] = ci_backscaled_2_2[1, :, 1]  # mac CI of shape 2 of ADE2
    final[:, 15] = labels_shape_r_2[:, 2]  # longevity of ADE2
    final[:, 16] = ci_backscaled_2_l[0, :, 2]  # min CI of longevity of ADE2
    final[:, 17] = ci_backscaled_2_l[1, :, 2]  # max CI of longevity of ADE2
    final[:, 18] = labels_shape_r_0[:, 0]  # shapes under no ADE
    final[:, 19] = ci_backscaled_0[0, :, 0]  # min CI no ADE
    final[:, 20] = ci_backscaled_0[1, :, 0]  # max CI no ADE
    return final, rmse_ade1, rmse_ade2_1, rmse_ade2_2

if __name__ == "__main__":
    #Load models
    wd_model = "/path_to_model_folder/"
    model_shape_1 = tf.keras.models.load_model(wd_model + "rnn_modelhistory_shape_1")
    model_single_shape = tf.keras.models.load_model(wd_model + "rnn_modelhistory_single_shape")
    model_two_shapes = tf.keras.models.load_model(wd_model + "rnn_modelhistory_two_shapes")
    model_clas = tf.keras.models.load_model(wd_model + "rnn_modelhistory_clas")

    #input file
    f = ("/path_to_input_file/Data_S1.txt")

    n_predictions = 100

    #Time bins determined by extinction patterns based on all neoselachii species
    # time_slice = [[145, 95.95], [95.95, 83.29], [84, 73.58], [73.58, 71.884], [66.782, 65.481], [65.481, 55.477], [55.477, 39.369], [39.369, 33.166], [33.166, 16.258], [16.258, 3.65], [3.65, 0.0117]]

    time_slice = [[145, 95.71], [95.71, 83.31], [83.31, 73.40], [73.40, 71.90], [66.60, 65.79],
                  [65.79, 55.89], [55.89, 38.58], [38.58, 33.57], [33.57, 3.65], [3.65, 0.0117]]  #  [71.90, 66.60] excluded due to low number of species

    # Create empty arrays to save the values
    data = np.zeros((len(time_slice), n_predictions, 21))
    rmse_1 = np.zeros(len(time_slice))
    rmse_2_1 = np.zeros(len(time_slice))
    rmse_2_2 = np.zeros(len(time_slice))

    for i in range(len(time_slice)):
        data[i, :, :], rmse_1[i,], rmse_2_1[i,], rmse_2_2[i,] = pipeline(f, n_predictions, time_slice[i])

    f1 = os.path.join("/path_to_adenn_output_folder/species_ADENN.npy")
    np.save(file=f1, arr=data)
