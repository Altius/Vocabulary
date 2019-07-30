# This script was used for the 2.5D NMF project. it uses data generated in the technical notebook "2-06-19 attempt to work with fimo - validation 7-30-19.ipynb"

import numpy as np
import pandas as pd
from OONMFhelpers import *

def acquire_data():
	finname = '20619_FIMOtable_clusters.txt'
	print('starting read at ', mytime(), flush=True)
	A = pd.read_table(finname)
	print('finished read at ', mytime(), flush=True)
	return A.values > 0.5 
	
A = acquire_data()
print('A shape',A.shape)
# should be 3.5e6 x 206

B = np.load('60518_NNDSVD_NC16/2018-06-08NC16_NNDSVD_Mixture.npy')
#should be 16 x 3.5e6

from datetime import date
from sklearn.decomposition import NMF, non_negative_factorization

today = str(date.today())

W, H, n_iter = non_negative_factorization(A.T, n_components=16, init='custom', random_state=3, update_H=False, H=B, verbose=1)

np.save(today+'MotifClusterPA_X_Components'+'NC16NNDSVD.npy', W)
