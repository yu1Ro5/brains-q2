import csv
import datetime
import os
import pickle
import statistics
import time
import warnings

import joblib
import lightgbm as lgb
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error as mse
from sklearn.model_selection import cross_validate
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR

import util


def lgbm(X,y,nf=5):
    dall=lgb.Dataset(data=X,label=y)
    best_params_650 = {
        'objective': 'regression',
        'metric': 'l2',
        'verbosity': -1,
        'boosting_type': 'gbdt',
        'force_col_wise': 'true',
        'num_threads': -1,
        'feature_pre_filter': False,
        'lambda_l1': 6.817063408573141e-07,
        'lambda_l2': 0.00012732475617557626,
        'num_leaves': 13,
        'feature_fraction': 0.5,
        'bagging_fraction': 0.9888992921439337,
        'bagging_freq': 1,
        'min_child_samples': 20
    }
    best_params={
        'objective': 'regression',
        'metric': 'l2',
        'verbosity': -1,
        'boosting_type': 'gbdt',
        'feature_pre_filter': False,
        'lambda_l1': 8.546942828085179e-06,
        'lambda_l2': 1.0468971074369936e-07,
        'num_leaves': 31,
        'feature_fraction': 0.4,
        'bagging_fraction': 1.0,
        'bagging_freq': 0,
        'min_child_samples': 20,
        'force_col_wise':'true'
    }
    lgb_cv = lgb.cv(
        params=best_params_650,
        train_set=dall,
        num_boost_round=1000,
        # folds=KFold(n_splits=5),
        nfold=nf,
        stratified=False,
        shuffle=True,
        metrics=["mean_squared_error"],
        early_stopping_rounds=100,
        seed=42,
        return_cvbooster=True
    )
    print("LGBM mean±std: ",lgb_cv['l2-mean'][-1],"±",lgb_cv['l2-stdv'][-1])
    joblib.dump(value=lgb_cv['cvbooster'].boosters,filename=open(file=os.path.dirname(__file__)+"/cv_estimators.joblib",mode="wb"),compress=9)
    joblib.dump(value=lgb_cv['cvbooster'].best_iteration,filename=open(file=os.path.dirname(__file__)+"/cv_best_iter.joblib",mode="wb"),compress=3)
    model_name = "LGBM"
    return lgb_cv['l2-mean'][-1],lgb_cv['l2-stdv'][-1],model_name

def xgbr(X,y):
    model_name = "XGBR"
    # xgb_params = {}
    XGBR = xgb.XGBRegressor(objective='reg:squarederror',verbosity=0,n_job=-1)
    # print(XGBR)
    cv_results = cross_validate(XGBR, X, y, cv=5, scoring='neg_mean_squared_error', return_estimator=True, n_jobs=-1)
    print("XGBR mean±std: ", -cv_results['test_score'].mean(),"±",statistics.pstdev(-cv_results['test_score']))
    joblib.dump(value=cv_results['estimator'],filename=open(file=os.path.dirname(__file__)+"/cv_estimators.joblib",mode="wb"),compress=9)
    return -cv_results['test_score'].mean(),statistics.pstdev(-cv_results['test_score']),model_name

def mlpr(X,y):
    model_name = "MLPR"
    # MLPR = MLPRegressor(hidden_layer_sizes=(200,200,200,100,),max_iter=1000)
    MLPR = MLPRegressor(activation='relu', alpha=0.0001, batch_size='auto', beta_1=0.9,
             beta_2=0.999, early_stopping=True, epsilon=1e-08,
             hidden_layer_sizes=(100, 100, 100,), learning_rate='constant',
             learning_rate_init=0.0016938817242064471, max_iter=10000,
             momentum=0.9, n_iter_no_change=10, nesterovs_momentum=True,
             power_t=0.5, random_state=42, shuffle=True, solver='adam',
             tol=0.0001, validation_fraction=0.1, verbose=False,
             warm_start=False)
    cv_results = cross_validate(MLPR, X, y, cv=5, scoring='neg_mean_squared_error', return_estimator=True, n_jobs=-1)
    # print(-cv_results['test_score'])
    print("MLPR mean±std: ", -cv_results['test_score'].mean(),"±",statistics.pstdev(-cv_results['test_score']))
    joblib.dump(value=cv_results['estimator'],filename=open(file=os.path.dirname(__file__)+"/cv_estimators.joblib",mode="wb"),compress=9)
    return -cv_results['test_score'].mean(),statistics.pstdev(-cv_results['test_score']),model_name

def ridge(X,y):
    model_name = "Ridge"
    RIDGE = Ridge()
    cv_results = cross_validate(RIDGE, X, y, cv=5, scoring='neg_mean_squared_error', return_estimator=True, n_jobs=-1)
    print("RIDGE mean±std: ", -cv_results['test_score'].mean(),"±",statistics.pstdev(-cv_results['test_score']))
    joblib.dump(value=cv_results['estimator'],filename=open(file=os.path.dirname(__file__)+"/cv_estimators.joblib",mode="wb"),compress=9)
    return -cv_results['test_score'].mean(),statistics.pstdev(-cv_results['test_score']),model_name

def svr(X,y,kernel='rbf'):
    _kernel = kernel
    model_name = "SVR " + _kernel
    svr_ = SVR(kernel=kernel, degree=3, gamma='auto')
    cv_results = cross_validate(svr_, X, y, cv=5, scoring='neg_mean_squared_error', return_estimator=True, n_jobs=-1)
    print("SVR mean±std: ", -cv_results['test_score'].mean(),"±",statistics.pstdev(-cv_results['test_score']))
    joblib.dump(value=cv_results['estimator'],filename=open(file=os.path.dirname(__file__)+"/cv_estimators.joblib",mode="wb"),compress=9)
    return -cv_results['test_score'].mean(),statistics.pstdev(-cv_results['test_score']),model_name

def log_train(now_time_str,X,prepared_time,cv_mean,cv_stdv,model_name):
    with open('log.csv', 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([now_time_str,X.shape[1],prepared_time,cv_mean,cv_stdv,model_name])

if __name__ == "__main__":
    warnings.simplefilter('ignore')
    start = time.time()
    df=pd.read_csv("datasets/dataset.csv")
    X,y=util.prepare_data(df)
    prepared_time = time.time() - start
    print ("prepared_time:{0}".format(prepared_time) + "[sec]")
    # print(X.shape,y.shape)
    # X=X.T.drop_duplicates().T
    print(X.shape,y.shape)
    cv_mean,cv_stdv,model_name = ridge(X,y)
    # cv_mean,cv_stdv,model_name = svr(X,y)
    cv_mean,cv_stdv,model_name = mlpr(X,y)
    # cv_mean,cv_stdv,model_name = lgbm(X,y)
    # cv_mean,cv_stdv,model_name = xgbr(X,y)
    # cv_mean,cv_stdv,model_name = ridge(X,y)
    now_time_str = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
    with open('log.csv', 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([now_time_str,X.shape[1],prepared_time,cv_mean,cv_stdv,model_name])
    li_model = [
        # svr(X,y),
        # mlpr(X,y),
        # lgbm(X,y),
        # lgbm(X,y,5,10),
        # lgbm(X,y,nf=5,md=7,nl=70),
        # ridge(X,y),
        # ori_cv(X,y,SVR()),
        # ori_cv(X,y,LGBMRegressor(),0),
        # ori_cv(X,y,LGBMRegressor(max_depth=7,num_leaves=30)),
        # ori_cv(X,y,LGBMRegressor(max_depth=9,num_leaves=70)),
        # ori_cv(X,y,LGBMRegressor(max_depth=14,num_leaves=160)),
        # ori_cv(X,y,MLPRegressor(max_iter=1000),1),
        # lgbm(X,y,nf=10,md=7,nl=30,save=True),
        # stack(X,y,1),
    ]

    # now_time_str = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
    # for i in li_model:
    #     log_train(now_time_str,X,prepared_time,i[0],i[1],i[2])
    #     if i == li_model[-1]:
    #         zip_filename = '.\\submission\\' + now_time_str + '_submit_'+str(X.shape[1])+'_'+i[2]+'.zip'
    #         print(zip_filename)
    #         with zipfile.ZipFile(zip_filename, 'w', compression=zipfile.ZIP_DEFLATED) as zip:
    #             add_all(zip, glob.glob("src/**", recursive=True))
    #             add_all(zip, glob.glob("env.yaml"))

    print("Good Job!")
