#%%
import pandas as pd
import scipy.io
import numpy as np
#import tensorflow as tftest4.py
import NeuralNetworkUtils as NN3

#path_db = '~/Users/jeanne.mermet-guyennet/Desktop/DTU/thesis/DB_case6/DB_final.mat'
#path_db = 'case14db.mat' # classification/classDetails
#DB = scipy.io.loadmat(path_db)
IN = np.loadtxt("./Last_dataset_from_Timon_18_Nov_2020/NN_input.csv", delimiter=",")#DB['classification']

OUT = np.loadtxt("./Last_dataset_from_Timon_18_Nov_2020/NN_output_mod.csv", delimiter=",")#DB['classification']
OUT = np.array([OUT])

DATA = np.concatenate((IN,OUT.T),axis=1)
#%% Build model
#                   M   H  H   H    N   eta,  phi, alpha
nn3 = NN3.Security([9, 50, 50, 50,  2], 0.001, 1, 0.001)

#%% Prepare data for testing
X_train, X_test, y_train, y_test, G_train, G_test = nn3.train_test_DB(DATA, test_size=0.2, save_splits=True)

#%% Prune model
crossval_3layer = nn3.train_and_prune(X_train, y_train, X_test, y_test, batch_size=500, nb_epochs=5000, disp=True, sparsity=.3, plot=True)

# Writing in files
path = './Last_dataset_from_Timon_18_Nov_2020/output/'
W1, b1, W2, b2, W3, b3,  W4, b4 = nn3.save_weights_bias(path)

