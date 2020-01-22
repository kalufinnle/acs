# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 17:46:58 2020

@author: kevin
"""
from sklearn.linear_model import Ridge
import numpy as np
from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_val_score

def load_peudo_data(n_samples, n_features):
    rng = np.random.RandomState(56)
    all_X = rng.randn(n_samples, n_features)
    all_y = rng.randn(n_samples)
    return all_X,all_y
    
def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * (np.sum(errors)/(errors.shape[0]))   
    accuracy = 100 - mape
    print('Model Performance')
    print('Test Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Test Accuracy = {:0.4f}%.'.format(accuracy))
    target_names = ['class 0', 'class 1']
    print(classification_report(test_labels, model.predict(test_features), target_names=target_names,digits=3))
    #confusion matrix and preparefor plot
    y_test = test_labels
    y_pred = model.predict(test_features)
    cnf_matrix = confusion_matrix(y_test, y_pred)
    print(str(cnf_matrix))
    return

def baseline_ridge(all_X,all_y):
    n_alphas = 200
    param_rdg = {'alpha':(np.logspace(-10, -2, n_alphas)).tolist()}
    ridge = linear_model.Ridge(fit_intercept=False)
    clf_s = GridSearchCV(ridge, param_grid = param_rdg, cv = 5)
    X_tr, X_ts, y_tr, y_ts = train_test_split(all_X, all_y, test_size = 0.2, random_state = 12, stratify = all_y)
    clf_s.fit(X_tr, y_tr)
    
    best_para = clf_s.best_params_
    best_model = clf_s.best_estimator_
    cv_for_hyperparatuning=clf_s.cv_results_
    ts_accuracy = evaluate(best_model, X_ts, y_ts)
    #run CV 100 times to smooth out the result//since reshuffle the data sometimes boost up the accuracy a lot   
    #cv = StratifiedShuffleSplit(n_splits=10, test_size=0.2, random_state=0)
    all_cvs=cross_val_score(best_model, X_tr, y_tr, cv = 5)  
    performance_cv_outcome = np.mean(all_cvs)
    print("CV Performance for Tuning =: %0.3f (+/- %0.3f)" % (all_cvs.mean(), all_cvs.std()))
    return cv_for_hyperparamtuning,best_para,best_model,ts_accuracy,performance_cv_outcome
    
def baseline_svc(all_X,all_y):
    param_svc = {'kernel':('linear','rbf','poly','sigmoid'), 'C':(np.linspace(0.1, 1, 100)).tolist()}
    svc = svm.SVC(gamma=0.001)
    clf_s = GridSearchCV(svc, param_grid = param_svc, cv = 5)
    X_tr, X_ts, y_tr, y_ts = train_test_split(all_X, all_y, test_size = 0.2, random_state = 12, stratify = all_y)
    #fit train data. Grid search CV is on train data //note that we have: train data/validation data/test data
    clf_s.fit(X_tr, y_tr) 
    
    best_para = clf_s.best_params_
    best_model = clf_s.best_estimator_
    cv_for_hyperparatuning=clf_s.cv_results_
    ts_accuracy = evaluate(best_model, X_ts, y_ts)
    #run CV 100 times to smooth out the result//since reshuffle the data sometimes boost up the accuracy a lot   
    #cv = StratifiedShuffleSplit(n_splits=10, test_size=0.2, random_state=0)
    all_cvs=cross_val_score(best_model, X_tr, y_tr, cv = 5)  
    performance_cv_outcome = np.mean(all_cvs)
    print("CV Performance for Tuning =: %0.3f (+/- %0.3f)" % (all_cvs.mean(), all_cvs.std()))
    return cv_for_hyperparamtuning,best_para,best_model,ts_accuracy,performance_cv_outcome