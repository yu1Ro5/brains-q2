import os
import pickle
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, PandasTools
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.preprocessing import StandardScaler

import util

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    input_data = []
    for line in sys.stdin:
        input_data.append(line.strip().split(","))

    df = pd.DataFrame(data=input_data[1:], columns=input_data[0])
    df = df.replace("", None)
    X,_ = util.prepare_data(df,False)
    X = X.astype(float)

    if 'ROMol' in X.columns:
        X = X.drop("ROMol",axis=1)

    cv_estimators = joblib.load(filename=open(file=os.path.dirname(__file__)+"/cv_estimators.joblib",mode="rb"))
    cv_bestiter = joblib.load(filename=open(file=os.path.dirname(__file__)+"/cv_best_iter.joblib",mode="rb"))

    preds = np.zeros(X.shape[0])
    for i in range(len(cv_estimators)):
        try:
            pred = cv_estimators[i].predict(data=X, num_iteration=cv_bestiter)
        except:
            pred = cv_estimators[i].predict(X)
        preds += pred
    preds = preds / len(cv_estimators)

    for val in preds:
        print(val)
