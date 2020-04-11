# 11-09-18 script to run NMF with new pancreatic samples added in. 

import numpy as np
import pandas as pd
import OONMF
from OONMFhelpers import *

today  = get_today()

def acquire_data():
        #fetch the presence absence masterlist
        print('starting read at ', mytime(), flush=True)
        A = pd.read_table('matrix_bin_all_733samples_WM20180608.txt', header=None)
        print('finished read at ', mytime(), flush=True)
        
        
        finname = 'alex_LN_names.txt'
        f = open(finname, 'r')
        dars = f.readlines()
        LN_ar = []
        for line in dars:
            LN_ar.append(line.strip())
        f.close()
        # these are the text vectors Shane helped Sasha generate on 11-08-18. one DHS per line, with 1 and 0 if the DHS is in a particular sample
        for LN in LN_ar:
            LNvals = np.loadtxt('final_Alex_peaks/ML_'+LN+'_E80p_overlapBMvec.txt').astype(int)
            #add these vectors to masterlist
            A[LN]=LNvals
        
        return A.drop([A.columns[0]], axis=1).values.T
        
# main masterlst presence absence matrix + new pancreas samples    
A = acquire_data()

#original 733 sample NMF  results for DHS component vectors
B = np.load('2018-06-08NC16_NNDSVD_Mixture.npy')

from datetime import date
from sklearn.decomposition import NMF, non_negative_factorization

today = str(date.today())

print('starting NMF at ',mytime())

#this is a special call for NMF which does not update one of the two matrices. In our case, we do not want to update the DHS x component matrix. 

W, H, n_iter = non_negative_factorization(A, n_components=16, init='custom', random_state=3, update_H=False, H=B, verbose=1,max_iter=20000)

print('done with NMF at ',mytime())
np.save(today+'PancSpecial'+'NC16seed20.npy', W)
