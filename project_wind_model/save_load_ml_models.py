# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:51:31 2017

@author: kevin

Saving Machine Learning Trained Models
"""

import pickle
'''e.g.
    from sklearn.linear_model import LogisticRegression
    import pickle
    model = LogisticRegression()
    model.fit(xtrain, ytrain)
'''

#deserialise save
def save_sklear_mod(model, model_file_path):
    # save the model to disk
    pickle.dump(model, open(model_file_path, 'wb'))
   
def load_sklear_mod(model_file_path):
    model = pickle.load(open(model_file_path, 'rb'))
    return model





    
#serialise save
def save_sklear_mod_s(model, model_file_path):
    # save the model to disk
    from sklearn.externals import joblib
    joblib.dump(model, model_file_path)
    
def load_sklear_mod_s(model_file_path):
    from sklearn.externals import joblib
    model = joblib.load(open(model_file_path, 'rb'))
    return model





  
def save_keras_mod(model, json_file_path, str_h5_file):  #str_h5_file = "model.h5"
    from keras.models import Sequential
    from keras.layers import Dense
    # serialize to JSON
    model_json = model.to_json()
    with open(json_file_path, "w") as file:
        file.write(model_json)
    # serialize weights to HDF5
    model.save_weights(str_h5_file)
    
def load_keras_mod(model, json_file_path, str_h5_file): #str_h5_file = "model.h5"
    from keras.models import model_from_json
    # load json
    file = open(json_file_path, 'r')
    model_json = file.read()
    file.close()
    loaded_model = model_from_json(model_json)
    # load weights
    loaded_model.load_weights(str_h5_file)
    return loaded_model




  
""" e.g.
    import statsmodels.api as sm
    model = sm.tsa.ARIMA([1,5,9,12], order=(1, 0, 1))  
    model= model.fit()
"""
def save_stats_mod(model, myfile):
    model.save(myfile)

def load_stats_mod(my_file):    
    from statsmodels.tsa.arima_model import ARIMAResults
    loaded = ARIMAResults.load(my_file)
    return loaded
