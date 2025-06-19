"""
Project: Diversification Rates and ADE
Description:
    Data simulation
    Model Training
    Model Accuracy, Coverage, RMSE
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

np.set_printoptions(suppress=True, precision=3)  # prints floats, no scientific notation
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

def print_update(s):
    sys.stdout.write('\r')
    sys.stdout.write(s)
    sys.stdout.flush()

def data_simulator_mixture(q=None,  # if None: rnd drawn
                           N=1000,
                           min_data_size=50,
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
                           mixture = False):
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
                x = np.random.randint(0,3)
                if x == 0: # both shapes == 1
                    rnd_shape = np.ones((1))
                    rnd_shape = np.tile(rnd_shape, (1, 2))
                    rnd_shapes[i] = rnd_shape
                elif x == 1: # both shapes identical
                    rnd_shape = np.exp(np.random.uniform(np.log(1 / magnitude), np.log(magnitude)))
                    rnd_shape= np.tile(rnd_shape, (1, 2))
                    rnd_shape = np.sort(rnd_shape, 1)
                    rnd_shapes[i] = rnd_shape
                elif x == 2:# two different shapes
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

    for i in np.arange(N):
        # loop over datasets
        EpochBins = np.linspace(0, 1000, 1 + int(1000 / EpochBinSize[i]))
        if i % 10 == 0:
            print_update("Simulating data...%s/%s" % (i, N))

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

        nF = np.random.poisson(q_rates[i] * rW)

        # actual occs| hist of how many sp extend over 1, 2, ... time bins
        rnd_spec_time = np.random.uniform(0, 10, n_species[i])
        hOccs = np.zeros(maxNbinsPerSpecies)
        brL = 0
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
            q_estimate = np.random.uniform(0.25, 1)
            sys.exit("???")
        else:
            q_estimate = np.sum(nF[nF > 1] - 1) / brL
        if verbose:
            print(q_rates[i], q_estimate)

        # BUILDING FEATURE SET
        # fraction of species with occurrences in 1, 2, 3, ... N time bins
        hOccs = hOccs / np.sum(hOccs)

        nF[nF > 1000] = 1000
        nFsampled = nF[nF > 0]
        if maxNbinsPerSpecies > 0:
            # append fraction of species with 1, 2, 3 ... N occurrences
            h_temp = np.append(np.histogram(nFsampled, bins=Hbins, density=True)[0], hOccs)
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
                "sample_species": len(nFsampled)
                }
    else:
        info = None
    return features, labels, info

def build_nn(dense_nodes=None,
             dense_act_f='relu',
             output_nodes=2,  # shape, mean longevity
             output_act_f='softplus',  # softplus/linear
             loss_f='mse',  # mean squared error
             dropout_rate=0,
             verbose=1):
    if dense_nodes is None:
        dense_nodes = [32]

    input_layer = [tf.keras.layers.Flatten(input_shape=[features.shape[1]])]
    model = keras.Sequential(input_layer)

    for i in range(len(dense_nodes)):
        model.add(layers.Dense(dense_nodes[i],
                               activation=dense_act_f))

    if dropout_rate:
        model.add(layers.Dropout(dropout_rate))

    model.add(layers.Dense(output_nodes,
                           activation=output_act_f))

    if verbose:
        print(model.summary())

    model.compile(loss=loss_f,
                  optimizer="adam",
                  metrics=[loss_f])
    return model


def fit_rnn(Xt, Yt, model,
            Xv=None, Yv=None,
            criterion="val_loss",
            patience=10,
            verbose=1,
            batch_size=1000,
            max_epochs=1000,
            validation_split=0.2,
            ):
    early_stop = keras.callbacks.EarlyStopping(monitor=criterion,
                                               patience=patience,
                                               restore_best_weights=True)

    if Xv is not None:
        validation_data = (Xv, Yv)
        validation_split = None

    history = model.fit(Xt, Yt,
                        epochs=max_epochs,
                        validation_split=validation_split,
                        verbose=verbose,
                        callbacks=[early_stop],
                        batch_size=batch_size,
                        validation_data=validation_data
                        )
    return history


def save_nn_model(wd, history, model, filename=""):
    # save training history
    with open(os.path.join(wd, "rnn_history" + filename + ".pkl"), 'wb') as output:  # Overwrites any existing file.
        pkl.dump(history.history, output, pkl.HIGHEST_PROTOCOL)
    # save model
    tf.keras.models.save_model(model, os.path.join(wd, 'rnn_model' + filename))


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

    def transform_labels(self, shapes, longevities):
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

        return tr_shapes, tr_longevities

    def back_transform_labels(self, tr_shapes, tr_longevities):
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

        return shapes, longevities


def plot_mixed_ade_predictions(true_lab, est_lab, alpha=0.5):
    # plot shapes
    fig = plt.figure(figsize=(14, 7))
    fig.add_subplot(1, 3, 1)
    plt.scatter(true_lab[:, 0], est_lab[:, 0], alpha=alpha)
    plt.scatter(true_lab[:, 1], est_lab[:, 1], alpha=alpha)
    plt.axline((0, 0), (1, 1), linewidth=2, linestyle='dashed', color="k")
    plt.xlabel('True shape')
    plt.ylabel('Estimated shape')

    # plot longevities
    fig.add_subplot(1, 3, 2)
    plt.scatter(true_lab[:, 2], est_lab[:, 2], alpha=alpha, )
    if true_lab.shape[1] == 4:
        plt.scatter(true_lab[:, 3], est_lab[:, 3], alpha=alpha, )
    plt.axline((0, 0), (1, 1), linewidth=2, linestyle='dashed', color="k")
    plt.xlabel('True longevity')
    plt.ylabel('Estimated longevity')

    # plot shape ratios
    fig.add_subplot(1, 3, 3)
    plt.scatter(true_lab[:, 1] / true_lab[:, 0], est_lab[:, 1] / est_lab[:, 0], alpha=alpha)
    plt.axline((0, 0), (1, 1), linewidth=2, linestyle='dashed', color="k")
    plt.xlabel('True shape ratio')
    plt.ylabel('Estimated shape ratio')

    plt.show()


def plot_ade_predictions(true_lab, est_lab, alpha=0.5):
    # plot shapes
    fig = plt.figure(figsize=(14, 7))
    fig.add_subplot(1, 2, 1)
    plt.scatter(true_lab[:, 0], est_lab[:, 0], alpha=alpha)
    plt.axline((0, 0), (1, 1), linewidth=2, linestyle='dashed', color="k")
    plt.xlabel('True shape')
    plt.ylabel('Estimated shape')

    # plot longevities
    fig.add_subplot(1, 2, 2)
    plt.scatter(true_lab[:, 1], est_lab[:, 1], alpha=alpha)
    plt.axline((0, 0), (1, 1), linewidth=2, linestyle='dashed', color="k")
    plt.xlabel('True longevity')
    plt.ylabel('Estimated longevity')

    plt.show()

def labels_to_class(labels_original, rescaler = None):
    class_labels = np.zeros((labels_original.shape[0], 3))
    labels = labels_original + 0
    if rescaler is not None:
        if rescaler.log_shapes:
            labels[:, 0] = np.exp(labels[:, 0])
            labels[:, 1] = np.exp(labels[:, 1])

    class_labels[np.where((labels[:, 0] == 1) & (labels[:, 1] == 1)), 0] = 1 #no ADE
    class_labels[np.where((labels[:, 0] == labels[:, 1]) & (labels[:, 0] != 1)), 1] = 1 #single ADE
    class_labels[np.where((np.sum(class_labels, 1) == 0)), 2] = 1 # double ADE

    return class_labels

def calcCI(predicted_labels, errors, mse = False):
    if mse == True:
        errors = np.sqrt(errors)
        CI1 = (predicted_labels - 1.96 * errors)
        CI2 = (predicted_labels + 1.96 * errors)
    else: #i.e. if errors are rmse
        CI1 = (predicted_labels - 1.96 * errors)
        CI2 = (predicted_labels + 1.96 * errors)
    return np.array([CI1, CI2])

if __name__ == '__main__':
    # simulate training set for ADE1 model
    # features, labels = data_simulator(N=10000, gamma_model=False)

    features, labels, info = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                    N=10000,
                                                    min_data_size=100,
                                                    max_data_size=1000,
                                                    fixed_shape=None,  # if not None -> array of 2 values
                                                    fixed_scale=None,  # if not None -> array of 2 values
                                                    magnitude=4,
                                                    # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                    min_longevity=2,
                                                    # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                    max_longevity=30,
                                                    # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                    n_Hbins=100,  # bins of histogram occs-per-species
                                                    maxNbinsPerSpecies=20,  # bins of histogram timebins-per-species
                                                    gamma_model=False,
                                                    alpha=1,
                                                    # rate heterogeneity across species (smaller alpha = greater variation)
                                                    verbose=0,
                                                    plot=False,
                                                    include_info=True,
                                                    single_longevity=True,
                                                    single_shape=True)

    # save to files
    wd = "/path_to_model_folder/"
    f1 = os.path.join(wd, "sim_features_mixture.npy")
    f2 = os.path.join(wd, "sim_labels.npy")
    np.save(file=f1, arr=features)
    np.save(file=f2, arr=labels)

    # features, labels = np.load(f1), np.load(f2)

    # data rescaler
    rescaler = AdeNNrescaler(log_shapes=True,
                             log_longevities=True,
                             shapes_mul=1,
                             longevities_mul=1,
                             shape_ratio=False)

    shapes_r, long_r = rescaler.transform_labels(labels[:, :2], labels[:, 2:])
    labels_r = np.hstack((shapes_r, long_r))

    features_train, features_val, labels_train, labels_val = train_test_split(features, labels_r,
                                                                              test_size=0.2, random_state=8)

    # build NN model
    model_ADE1 = build_nn(dense_nodes=[64, 8],
                     output_nodes=labels.shape[1],
                     output_act_f="linear",
                     dropout_rate=0.05)

    # train NN
    history = fit_rnn(features_train, labels_train, model, features_val, labels_val,
                      batch_size=100, max_epochs=3000)

    n_predictions = 1000
    labels_val_predict = np.array([model_ADE1(features_val, training=True) for _ in range(n_predictions)])
    labels_val_predict_mean = np.mean(labels_val_predict, axis=0)

    f3 = os.path.join(wd, "labels_predicted_val.npy")
    np.save(file=f3, arr=labels_val_predict)

    # plot predictions (transformed values)
    plot_mixed_ade_predictions(labels_val, labels_val_predict_mean)

    # plot predicitons (original values)
    pred_shapes, pred_long = rescaler.back_transform_labels(labels_val_predict_mean[:, :2],
                                                            labels_val_predict_mean[:, 2:])
    pred_labels = np.hstack((pred_shapes, pred_long))

    val_shapes, val_long = rescaler.back_transform_labels(labels_val[:, :2], labels_val[:, 2:])
    val_labels = np.hstack((val_shapes, val_long))

    plot_mixed_ade_predictions(val_labels, pred_labels)

    # sims with shape ratio > 2
    indx = np.where(val_labels[:, 1] / val_labels[:, 0] > 2)
    plot_mixed_ade_predictions(val_labels[indx], pred_labels[indx], alpha=0.5)

    # sims with shape n.2 > 2
    indx = np.where(val_labels[:, 1] > 2)
    plot_mixed_ade_predictions(val_labels[indx], pred_labels[indx], alpha=0.5)

    # sims with shape n.1 > 2
    indx = np.where(val_labels[:, 0] > 2)
    plot_mixed_ade_predictions(val_labels[indx], pred_labels[indx], alpha=0.5)

    #calculate and save model rmse
    rmse = np.sqrt(np.mean((val_labels - pred_labels) ** 2, axis=0))
    f4 = os.path.join(wd, "rmse.npy")
    np.save(file=f4, arr=rmse)

    # save model
    save_nn_model(wd, history, model, filename="history_shape_1")


   # simulate test set for ADE1 model
    features_test, labels_test, info_test = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                                   N=1000,
                                                                   min_data_size=100,
                                                                   max_data_size=1000,
                                                                   fixed_shape=None,  # if not None -> array of 2 values
                                                                   fixed_scale=None,  # if not None -> array of 2 values
                                                                   magnitude=4,
                                                                   # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                                   min_longevity=2,
                                                                   # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                                   max_longevity=30,
                                                                   # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                                   n_Hbins=100,  # bins of histogram occs-per-species
                                                                   maxNbinsPerSpecies=20,
                                                                   # bins of histogram timebins-per-species
                                                                   gamma_model=False,
                                                                   alpha=1,
                                                                   # rate heterogeneity across species (smaller alpha = greater variation)
                                                                   verbose=0,
                                                                   plot=False,
                                                                   include_info=True,
                                                                   single_longevity=True,
                                                                   single_shape=True,
                                                                   )

    #save test datasets
    f5 = os.path.join(wd, "sim_features_test_mixture.npy")
    f6 = os.path.join(wd, "sim_labels_test.npy")
    np.save(file=f5, arr=features_test)
    np.save(file=f6, arr=labels_test)

    shapes_test_r, long_test_r = rescaler.transform_labels(labels_test[:, :2], labels_test[:, 2:])
    labels_test_r = np.hstack((shapes_test_r, long_test_r))

    n_predictions = 1000  #
    pred_test = np.array([model_ADE1(features_test, training=True) for _ in range(n_predictions)])


    pred_test_mean = np.mean(pred_test, axis=0)

    f7 = os.path.join(wd, "labels_predicted_transformed.npy")
    np.save(file=f7, arr=pred_test_mean)

    plot_mixed_ade_predictions(labels_test_r, pred_test_mean)

    pred_test_shapes, pred_test_long = rescaler.back_transform_labels(pred_test_mean[:, :2], pred_test_mean[:, 2:])
    pred_test_labels = np.hstack((pred_test_shapes, pred_test_long))

    f8 = os.path.join(wd, "labels_predicted_original.npy")
    np.save(file=f8, arr=pred_test_labels)

    plot_mixed_ade_predictions(labels_test, pred_test_labels)

    #Calculate confidence intervals
    confidence_intervals_shape_1 = np.zeros((labels_test.shape[0], 2))
    confidence_intervals_longevity = np.zeros((labels_test.shape[0], 2))
    for i in range(labels_test.shape[0]):
        confidence_intervals_shape_1[i, ] = calcCI(pred_test_labels[i, 0], errors=rmse[0]) #skipping for second shape as it is identical to first shape
        confidence_intervals_longevity[i, ] = calcCI(pred_test_labels[i, 2], errors=rmse[2])

    #Calculate coverage

    hits_shape_1 = np.asarray(np.where((labels_test[:, 0] >= confidence_intervals_shape_1[:, 0]) & (labels_test[:, 0] <= confidence_intervals_shape_1[:, 1]))).shape[1]
    hits_longevity = np.asarray(np.where((labels_test[:, 2] >= confidence_intervals_longevity[:, 0]) & (labels_test[:, 2] <= confidence_intervals_longevity[:, 1]))).shape[1]

    coverage_shape_1 = (hits_shape_1/labels_test.shape[0])*100
    coverage_longevity = (hits_longevity / labels_test.shape[0]) * 100




    #Simulate training data for ADE2 model

    features, labels, info = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                    N=10000,
                                                    min_data_size=100,
                                                    max_data_size=1000,
                                                    fixed_shape=None,  # if not None -> array of 2 values
                                                    fixed_scale=None,  # if not None -> array of 2 values
                                                    magnitude=4,
                                                    # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                    min_longevity=2,
                                                    # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                    max_longevity=30,
                                                    # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                    n_Hbins=100,  # bins of histogram occs-per-species
                                                    maxNbinsPerSpecies=20,  # bins of histogram timebins-per-species
                                                    gamma_model=False,
                                                    alpha=1,
                                                    # rate heterogeneity across species (smaller alpha = greater variation)
                                                    verbose=0,
                                                    plot=False,
                                                    include_info=True,
                                                    single_longevity=True,
                                                    single_shape=False,
                                                    two_shapes= True)

    # save to files
    wd = "/path_to_model_folder/"
    f1 = os.path.join(wd, "sim_features_mixture.npy")
    f2 = os.path.join(wd, "sim_labels.npy")
    np.save(file=f1, arr=features)
    np.save(file=f2, arr=labels)

    # features, labels = np.load(f1), np.load(f2)

    # data rescaler
    rescaler = AdeNNrescaler(log_shapes=True,
                             log_longevities=True,
                             shapes_mul=1,
                             longevities_mul=1,
                             shape_ratio=False)

    shapes_r, long_r = rescaler.transform_labels(labels[:, :2], labels[:, 2:])
    labels_r = np.hstack((shapes_r, long_r))

    features_train, features_val, labels_train, labels_val = train_test_split(features, labels_r,
                                                                              test_size=0.2, random_state=8)

    # build NN model
    model_ADE2 = build_nn(dense_nodes=[64, 8],
                          output_nodes=labels.shape[1],
                          output_act_f="linear",
                          dropout_rate=0.05)

    # train NN
    history = fit_rnn(features_train, labels_train, model, features_val, labels_val,
                      batch_size=100, max_epochs=3000)

    n_predictions = 1000
    labels_val_predict = np.array([model_ADE2(features_val, training=True) for _ in range(n_predictions)])
    labels_val_predict_mean = np.mean(labels_val_predict, axis=0)

    f3 = os.path.join(wd, "labels_predicted_val.npy")
    np.save(file=f3, arr=labels_val_predict)

    # plot predictions (transformed values)
    plot_mixed_ade_predictions(labels_val, labels_val_predict_mean)

    # plot predicitons (original values)
    pred_shapes, pred_long = rescaler.back_transform_labels(labels_val_predict_mean[:, :2],
                                                            labels_val_predict_mean[:, 2:])
    pred_labels = np.hstack((pred_shapes, pred_long))

    val_shapes, val_long = rescaler.back_transform_labels(labels_val[:, :2], labels_val[:, 2:])
    val_labels = np.hstack((val_shapes, val_long))

    plot_mixed_ade_predictions(val_labels, pred_labels)

    # sims with shape ratio > 2
    indx = np.where(val_labels[:, 1] / val_labels[:, 0] > 2)
    plot_mixed_ade_predictions(val_labels[indx], pred_labels[indx], alpha=0.5)

    # sims with shape n.2 > 2
    indx = np.where(val_labels[:, 1] > 2)
    plot_mixed_ade_predictions(val_labels[indx], pred_labels[indx], alpha=0.5)

    # sims with shape n.1 > 2
    indx = np.where(val_labels[:, 0] > 2)
    plot_mixed_ade_predictions(val_labels[indx], pred_labels[indx], alpha=0.5)

    # calculate and save model rmse
    rmse = np.sqrt(np.mean((val_labels - pred_labels) ** 2, axis=0))
    f4 = os.path.join(wd, "rmse.npy")
    np.save(file=f4, arr=rmse)

    # save model
    save_nn_model(wd, history, model, filename="history_shape_1")

    # simulate test set for ADE2 model
    features_test, labels_test, info_test = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                                   N=1000,
                                                                   min_data_size=100,
                                                                   max_data_size=1000,
                                                                   fixed_shape=None,  # if not None -> array of 2 values
                                                                   fixed_scale=None,  # if not None -> array of 2 values
                                                                   magnitude=4,
                                                                   # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                                   min_longevity=2,
                                                                   # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                                   max_longevity=30,
                                                                   # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                                   n_Hbins=100,  # bins of histogram occs-per-species
                                                                   maxNbinsPerSpecies=20,
                                                                   # bins of histogram timebins-per-species
                                                                   gamma_model=False,
                                                                   alpha=1,
                                                                   # rate heterogeneity across species (smaller alpha = greater variation)
                                                                   verbose=0,
                                                                   plot=False,
                                                                   include_info=True,
                                                                   single_longevity=True,
                                                                   single_shape=False,
                                                                   two_shapes= True
                                                                   )

    # save test datasets
    f5 = os.path.join(wd, "sim_features_test_mixture.npy")
    f6 = os.path.join(wd, "sim_labels_test.npy")
    np.save(file=f5, arr=features_test)
    np.save(file=f6, arr=labels_test)

    shapes_test_r, long_test_r = rescaler.transform_labels(labels_test[:, :2], labels_test[:, 2:])
    labels_test_r = np.hstack((shapes_test_r, long_test_r))

    n_predictions = 1000  #
    pred_test = np.array([model_ADE2(features_test, training=True) for _ in range(n_predictions)])

    pred_test_mean = np.mean(pred_test, axis=0)

    f7 = os.path.join(wd, "labels_predicted_transformed.npy")
    np.save(file=f7, arr=pred_test_mean)

    plot_mixed_ade_predictions(labels_test_r, pred_test_mean)

    pred_test_shapes, pred_test_long = rescaler.back_transform_labels(pred_test_mean[:, :2], pred_test_mean[:, 2:])
    pred_test_labels = np.hstack((pred_test_shapes, pred_test_long))

    f8 = os.path.join(wd, "labels_predicted_original.npy")
    np.save(file=f8, arr=pred_test_labels)

    plot_mixed_ade_predictions(labels_test, pred_test_labels)

    # Calculate confidence intervals
    confidence_intervals_shape_1 = np.zeros((labels_test.shape[0], 2))
    confidence_intervals_shape_2 = np.zeros((labels_test.shape[0], 2))
    confidence_intervals_longevity = np.zeros((labels_test.shape[0], 2))
    for i in range(labels_test.shape[0]):
        confidence_intervals_shape_1[i,] = calcCI(pred_test_labels[i, 0], errors=rmse[
            0])
        confidence_intervals_shape_2[i,] = calcCI(pred_test_labels[i, 1], errors=rmse[
            1])
        confidence_intervals_longevity[i,] = calcCI(pred_test_labels[i, 2], errors=rmse[2])

    # Calculate coverage

    hits_shape_1 = np.asarray(np.where((labels_test[:, 0] >= confidence_intervals_shape_1[:, 0]) & (
                labels_test[:, 0] <= confidence_intervals_shape_1[:, 1]))).shape[1]
    hits_shape_2 = np.asarray(np.where((labels_test[:, 1] >= confidence_intervals_shape_2[:, 0]) & (
                labels_test[:, 1] <= confidence_intervals_shape_2[:, 1]))).shape[1]
    hits_longevity = np.asarray(np.where((labels_test[:, 2] >= confidence_intervals_longevity[:, 0]) & (
                labels_test[:, 2] <= confidence_intervals_longevity[:, 1]))).shape[1]

    coverage_shape_1 = (hits_shape_1 / labels_test.shape[0]) * 100
    coverage_shape_2 = (hits_shape_2 / labels_test.shape[0]) * 100
    coverage_longevity = (hits_longevity / labels_test.shape[0]) * 100

    # Simulate training data for classifier model

    features, labels, info = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                    N=10000,
                                                    min_data_size=100,
                                                    max_data_size=1000,
                                                    fixed_shape=None,  # if not None -> array of 2 values
                                                    fixed_scale=None,  # if not None -> array of 2 values
                                                    magnitude=4,
                                                    # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                    min_longevity=2,
                                                    # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                    max_longevity=30,
                                                    # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                    n_Hbins=100,  # bins of histogram occs-per-species
                                                    maxNbinsPerSpecies=20,  # bins of histogram timebins-per-species
                                                    gamma_model=False,
                                                    alpha=1,
                                                    # rate heterogeneity across species (smaller alpha = greater variation)
                                                    verbose=0,
                                                    plot=False,
                                                    include_info=True,
                                                    single_longevity=True,
                                                    single_shape=False,
                                                    two_shapes=False,
                                                    mixture= True)

    # save to files
    wd = "/path_to_model_folder/single_shape/training/"
    f1 = os.path.join(wd, "sim_features_mixture.npy")
    f2 = os.path.join(wd, "sim_labels.npy")
    np.save(file=f1, arr=features)
    np.save(file=f2, arr=labels)

    # features, labels = np.load(f1), np.load(f2)

    model_clas = build_nn(dense_nodes=[128, 64, 64],
                          output_nodes=3,
                          output_act_f="softmax",  # classifier function
                          loss_f='categorical_crossentropy',
                          )

    shapes_r, long_r = rescaler.transform_labels(labels[:, :2], labels[:, 2:])
    labels_r = np.hstack((shapes_r, long_r))

    features_train, features_val, labels_train, labels_val = train_test_split(features, labels_r,
                                                                              test_size=0.2, random_state=8)

    history = fit_rnn(features_train, labels_to_class(labels_train, rescaler), model_clas, features_val,
                      labels_to_class(labels_val, rescaler),
                      batch_size=1000, max_epochs=3000)
    save_nn_model(wd, history, model_clas, filename="history_clas")
    # when fitting model_clas apply labels_to_class on both train and validation data

    # simulate test set for classifier model
    features_test, labels_test, info_test = data_simulator_mixture(q=None,  # if None: rnd drawn
                                                                   N=1000,
                                                                   min_data_size=100,
                                                                   max_data_size=1000,
                                                                   fixed_shape=None,  # if not None -> array of 2 values
                                                                   fixed_scale=None,  # if not None -> array of 2 values
                                                                   magnitude=4,
                                                                   # shape range: exp( U[log(1/magnitude), log(magnitude)] )
                                                                   min_longevity=2,
                                                                   # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                                   max_longevity=30,
                                                                   # ~ exp( U[log(min_longevity), log(max_longevity)] )
                                                                   n_Hbins=100,  # bins of histogram occs-per-species
                                                                   maxNbinsPerSpecies=20,
                                                                   # bins of histogram timebins-per-species
                                                                   gamma_model=False,
                                                                   alpha=1,
                                                                   # rate heterogeneity across species (smaller alpha = greater variation)
                                                                   verbose=0,
                                                                   plot=False,
                                                                   include_info=True,
                                                                   single_longevity=True,
                                                                   single_shape=False,
                                                                   two_shapes= False,
                                                                   mixture= True
                                                                   )
    wd = "/path_to_model_folder/single_shape/test/"
    f1 = os.path.join(wd, "sim_features_mixture_test.npy")
    f2 = os.path.join(wd, "sim_labels_test.npy")
    np.save(file=f1, arr=features)
    np.save(file=f2, arr=labels)

    rescaler = AdeNNrescaler(log_shapes=True,
                             log_longevities=True,
                             shapes_mul=1,
                             longevities_mul=1,
                             shape_ratio=False)

    pred_test = model_clas.predict(features_test)
    Truth = np.argmax(labels_to_class(labels_test), axis=1)
    # Truth = np.argmax(labels_to_class(labels_train, rescaler), axis = 1)
    Predict = np.argmax(pred_test, axis=1)

    conf_mat = confusion_matrix(Truth, Predict)
    sns.heatmap(conf_mat, cmap="rocket_r")

    #Calculate accuracy of the classifier model
    sums = np.sum(conf_mat, axis=1) #true numbers of classes
    correct = np.array([conf_mat[0, 0], conf_mat[1, 1], conf_mat[2, 2]]) #correctly classified classes

    total = (np.sum(correct) / np.sum(sums)) * 100
    individual_classes = (correct / sums) * 100

    #Calculate false positives (ADE0 gets rejected when true)

rescaler = AdeNNrescaler(log_shapes=True,
                         log_longevities=True,
                         shapes_mul=1,
                         longevities_mul=1,
                         shape_ratio=False)

features_test = np.load("/path_to_model_folder/sim_features_500.npy")
labels_test = np.load("/path_to_model_folder/sim_labels_500.npy")

pred_test = model_clas.predict(features_test)
Truth = np.argmax(labels_to_class(labels_test), axis=1)
# Truth = np.argmax(labels_to_class(labels_train, rescaler), axis = 1)
Predict = np.argmax(pred_test, axis=1)

indx = np.where(Predict == 1)

features_test_filtered = features_test[indx]

n_predictions = 1000  #

labels_filtered_predicted = np.array([model_single_shape(features_test_filtered) for _ in range(n_predictions)])
pred_test_mean = np.mean(labels_filtered_predicted, axis=0)

pred_test_shapes, pred_test_long = rescaler.back_transform_labels(pred_test_mean[:, :2], pred_test_mean[:, 2:])
pred_test_labels = np.hstack((pred_test_shapes, pred_test_long))

#get rmse for the appropriare sample size

def get_rmse(sample_size, model):
    if model == "ADE1":
        features, labels, info = data_simulator_mixture(q=1,  # if None: rnd drawn
                                                        N=1000,
                                                        min_data_size=sample_size - 1,
                                                        max_data_size=sample_size,
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
        rescaler = AdeNNrescaler(log_shapes=True,
                             log_longevities=True,
                             shapes_mul=1,
                             longevities_mul=1,
                             shape_ratio=False)
        shapes_r, long_r = rescaler.transform_labels(labels[:, :2], labels[:, 2:])
        labels_r = np.hstack((shapes_r, long_r))

        n_predictions = 1000  #
        pred = np.array([model_single_shape(features, training=True) for _ in range(n_predictions)])
        pred_mean = np.mean(pred, axis=0)

        rmse = np.sqrt(np.mean((labels_r - pred_mean) ** 2, axis=0))

    elif model == "ADE2":
        features, labels, info = data_simulator_mixture(q=1,  # if None: rnd drawn
                                                        N=1000,
                                                        min_data_size=sample_size - 1,
                                                        max_data_size=sample_size,
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
        rescaler = AdeNNrescaler(log_shapes=True,
                                 log_longevities=True,
                                 shapes_mul=1,
                                 longevities_mul=1,
                                 shape_ratio=False)
        shapes_r, long_r = rescaler.transform_labels(labels[:, :2], labels[:, 2:])
        labels_r = np.hstack((shapes_r, long_r))

        n_predictions = 1000  #
        pred = np.array([model_two_shapes(features, training=True) for _ in range(n_predictions)])
        pred_mean = np.mean(pred, axis=0)

        rmse = np.sqrt(np.mean((labels_r - pred_mean) ** 2, axis=0))

    return rmse

rmse_ADE1 = get_rmse(500, "ADE1")
rmse_ADE2 = get_rmse(500, "ADE2")

wd = "/path_to_model_folder/"
f1 = os.path.join(wd, "rmse_ADE1_500.npy")
np.save(file=f1, arr=rmse_ADE1)

rmse_ADE1 = np.load("/path_to_model_folder/rmse_ADE1_500.npy")
rmse_ADE2 = np.load("/path_to_model_folder/rmse_ADE2_500.npy")

#for ADE1

CIs = calcCI(pred_test_labels, errors=rmse_ADE1)

counter_ADE1 = 0

for i in range(CIs.shape[1]):
    if CIs[0, i, 0] < 1 and CIs[1, i, 0] > 1:
        counter_ADE1 +=1

false_positives_ADE1 = CIs.shape[1] - counter_ADE1

#for ADE2

CIs_1 = calcCI(pred_test_labels[:, 0], errors=rmse_ADE2[0])
CIs_2 = calcCI(pred_test_labels[:,1], errors=rmse_ADE2[1])

counter_ADE2 = 0

for i in range(CIs_1.shape[1]):
    if CIs_1[0, i] < 1 and CIs_1[1, i] > 1:
        if CIs_2[0, i] < 1 and CIs_2[1, i] > 1:
            counter_ADE2 += 1

false_positives_ADE2 = CIs_1.shape[1] - counter_ADE2
