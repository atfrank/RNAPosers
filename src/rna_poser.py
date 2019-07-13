from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import scipy.stats as stats
import numpy, string
import os
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--classifier", help="classifier")
parser.add_argument("--features", help="molecularize features: e.g: /home/itssahil/PROJECTS/Mol_Feturizer/B_LATEST_results/PREDICTORS_R/RF_6bfb.predictor.pkl")
parser.add_argument("--output", help="output scores")

a = parser.parse_args()

def main():
    MAXF, MINF = 3.4e38, 1.18e-38
    Xtest = numpy.loadtxt(a.features, dtype='float64')
    Xtest = Xtest[:, 1:]
    # workaround overflow
    Xtest = numpy.nan_to_num(Xtest)
    Xtest[Xtest < MINF] = MINF
    Xtest[Xtest > MAXF] = MAXF
    # predict
    pred_rf = joblib.load(a.classifier)
    rf = pred_rf.predict(Xtest)
    rf_p = pred_rf.predict_proba(Xtest)
    merge_rf = numpy.column_stack((rf, rf_p))
    numpy.savetxt(a.output, merge_rf, fmt='%f')
main()
