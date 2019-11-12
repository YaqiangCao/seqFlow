#!/usr/bin/env python2.7
#--coding:utf-8--
"""
Deep-learning neuron network model for richs prediction, richs were detected from Trac-looping data. 
The model is also used for feature importance detection.
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-02-09"
__modified__ = "2015-08-14"
__email__ = "caoyaqiang0410@gmail.com"

#systematic library
import glob, os, time
from datetime import datetime
from copy import deepcopy

#computating setting
import pandas as pd
import numpy as np
from tqdm import tqdm

#machine learning related
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.utils import class_weight

#deep-learning related
import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import load_model, Sequential, Model
from tensorflow.keras.callbacks import ModelCheckpoint, ReduceLROnPlateau, EarlyStopping
from tensorflow.keras.layers import Flatten, Dense, Dropout, BatchNormalization, Activation

#cLoops2 settings
from cLoops2.settings import *


#gloabl settings
class hyperparameters:
    num_classes = 3
    dim = 10
    batch_size = 2048
    learning_rate = 1e-2
    weight_decay = 0.0005
    train_steps = 1000  # will be changed accoridng to sample numbers
    vali_steps = 100
    test_steps = 100
    epochs = 100
    reduce_lr_patience = 10
    early_stop_patience = 20
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = '3'


PARA = hyperparameters()


def readData(f, ylabel="rich"):
    """
    Prepare the data for trainning.
    """
    mat = pd.read_csv(f, index_col=0, sep="\t")
    mat = mat.fillna(0.0)
    y = mat[ylabel]
    y = y.astype(int)
    del mat[ylabel]
    x = mat
    return x, y


def getY(y, classes):
    y = y.values
    ny = []
    for i in range(len(y)):
        n = [0] * classes
        n[y[i]] = 1
        ny.append(n)
    ny = np.array(ny)
    return ny


def getModel(input_shape,
             num_classes,
             checkpoint=None,
             lr=1e-3,
             reduce_lr_patience=5,
             early_stop_patience=10):
    if os.path.isfile(checkpoint):
        print("loading existing model")
        model = load_model(checkpoint)
    else:
        model = Sequential()
        model.add(Dense(64, input_dim=input_shape, activation="relu"))
        model.add(Dense(32, activation="relu"))
        model.add(Dense(16, activation="relu"))
        model.add(Dropout(0.1))
        model.add(Dense(num_classes, activation="softmax"))
        model.compile(
            loss="binary_crossentropy",
            optimizer=Adam(lr=lr),
            metrics=['accuracy'],
        )
    reduce_lr = ReduceLROnPlateau(monitor="val_loss",
                                  patience=reduce_lr_patience)
    early_stop = EarlyStopping(monitor='val_loss',
                               patience=early_stop_patience)
    callbacks = [reduce_lr, early_stop]
    if checkpoint is not None:
        cp = ModelCheckpoint(checkpoint,
                             monitor='val_loss',
                             verbose=1,
                             save_weights_only=False,
                             save_best_only=True,
                             mode='min')
    callbacks.append(cp)
    return model, callbacks


def getStat(model, x, y):
    yps = []
    yps = model.predict(x)
    yps = [np.argmax(yp) for yp in yps]
    acc = accuracy_score(y, yps)
    return acc


def train(trainF, valiF, checkpoint):
    x_train, y_train = readData(trainF)
    x_vali, y_vali = readData(valiF)
    PARA.dim = x_train.shape[1]
    PARA.num_classes = len(set(y_train))
    PARA.train_steps = len(y_train) / PARA.batch_size
    PARA.vali_steps = len(y_vali) / PARA.batch_size
    class_weights = class_weight.compute_class_weight("balanced",
                                                      np.unique(y_train),
                                                      y_train)
    y_train = getY(y_train, PARA.num_classes)
    y_vali = getY(y_vali, PARA.num_classes)
    model, callbacks = getModel(PARA.dim,
                                PARA.num_classes,
                                checkpoint=checkpoint,
                                reduce_lr_patience=PARA.reduce_lr_patience,
                                early_stop_patience=PARA.early_stop_patience)
    hist = model.fit(x=x_train,
                     y=y_train,
                     callbacks=callbacks,
                     epochs=PARA.epochs,
                     shuffle=True,
                     class_weight=class_weights,
                     validation_data=(x_vali, y_vali))
    K.clear_session()
    hist = pd.DataFrame(hist.history)
    hist.to_csv("trainningHistroy.txt", sep="\t", index_label="epoch")
    K.clear_session()


def getFeatureImportance(cp, f, pre="", repeat=10):
    """
    Get the feature importance by permutation input features
    """
    model = load_model(cp)
    x, y = readData(f)
    yps = model.predict(x)
    yps = [np.argmax(yp) for yp in yps]
    rawacc = accuracy_score(y, yps)  #raw acc used to compare decrease of acc
    ss = {}
    print("starting evulating feature importance")
    for c in tqdm(x.columns):
        rawc = deepcopy(x[c])
        ss[c] = {}
        for i in range(repeat):
            newc = deepcopy(x[c].values)
            np.random.shuffle(newc)
            newc = pd.Series(newc, index=rawc.index)
            x[c] = newc
            yps = model.predict(x)
            yps = [np.argmax(yp) for yp in yps]
            acc = accuracy_score(y, yps)
            ss[c][i] = rawacc - acc
        x[c] = rawc
    ss = pd.DataFrame(ss).T
    ss.to_csv("%s_rich_feature_importance.txt" % pre, sep="\t")
    ns = ss.mean(axis=1)
    ns = ns.sort_values(inplace=False, ascending=False)


def plotTrainning(trainF, valiF, testF, cp):
    x_train, y_train = readData(trainF)
    x_vali, y_vali = readData(valiF)
    x_test, y_test = readData(testF)
    model = load_model(cp)
    acc_train = getStat(model, x_train, y_train)
    acc_vali = getStat(model, x_vali, y_vali)
    acc_test = getStat(model, x_test, y_test)
    mat = pd.read_csv("trainningHistroy.txt", index_col=0, sep="\t")
    fig, ax = pylab.subplots()
    ax.plot(mat.index,
            mat["acc"],
            color=colors[0],
            linewidth=2,
            label="trainning:%.3f acc" % acc_train)
    ax.plot(mat.index,
            mat["val_acc"],
            color=colors[1],
            linewidth=2,
            label="validation:%.3f acc" % acc_vali)
    ax.set_xlabel("epoch")
    ax.set_ylabel("accuracy")
    ax.set_ylim([0.8, 1.0])
    ax.legend()
    ax.set_title("test set:%.3f acc" % acc_test)
    pylab.savefig("trainningHistroy.pdf")


def plotValiAUC(valiF, testF, cp):
    model = load_model(cp)
    x_vali, y_vali = readData(valiF)
    yps = model.predict(x_vali)
    yps = yps[:, 0]
    fpr, tpr, thresholds = roc_curve(y_vali, yps, pos_label=0)
    auc_vali = auc(fpr, tpr)
    x_test, y_test = readData(testF)
    yps = model.predict(x_test)
    yps = yps[:, 0]
    auc_test = auc(fpr, tpr)
    fig, ax = pylab.subplots(figsize=(2 * 0.8, 2.75 * 0.8))
    ax.bar([0, 1], [auc_vali, auc_test], color=[colors[0], colors[1]])
    ax.set_xticklabels(["vali", "test"])
    ax.set_ylabel("AUC")
    ax.set_title("AUC for RICH classificaiton")
    ax.set_ylim([0.9, 1.0])
    pylab.savefig("roc.pdf")


def plotFeatureImportance(f, pre):
    mat = pd.read_csv(f, index_col=0, sep="\t")
    ss = mat.mean(axis=1)
    ss = ss.sort_values(inplace=False, ascending=False)
    print(ss[:20])


def main():
    trainF = "../5.1.data/richTrain.txt"
    valiF = "../5.1.data/richVali.txt"
    testF = "../5.1.data/richTest.txt"
    cp = "richModel.h5"
    train(trainF, valiF, cp)
    getFeatureImportance(cp, valiF, "vali")
    getFeatureImportance(cp, testF, "test")
    plotTrainning(trainF, valiF, testF, cp)
    plotValiAUC(valiF, testF, cp)
    plotFeatureImportance("vali_rich_feature_importance.txt", "vali")
    plotFeatureImportance("test_rich_feature_importance.txt", "test")


if __name__ == "__main__":
    main()
