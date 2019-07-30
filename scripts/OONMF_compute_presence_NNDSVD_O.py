#libraries
import numpy as np
import pandas as pd
import OONMF
from OONMFhelpers import *

today  = get_today()

#presence absence matrix for 733 ENCODE samples formatted in a specific way
def acquire_data():
	print('starting read at ', mytime(), flush=True)
	A = pd.read_table('matrix_bin_all_733samples_WM20180608.txt', header=None)
	print('finished read at ', mytime(), flush=True)
	return A.drop([A.columns[0]], axis=1).values.T
	
A = acquire_data()

#number of components (k)
Nc = 16

#random seed (not very important for NNDSVD)
seed = 20

#instantiate object
a = OONMF.NMFobject(theNcomps=Nc)

#perform NMF with NNDSVD.  Requires about 50-100 GB memory. Should take 1 hour on Altius HPC cluster. 
a.performNMF(data=A, randomseed=seed, theinit='nndsvd')	

#write the output matrices into numpy binary
a.writeNMF(Basis_foutname= today+'NC'+str(Nc)+'_NNDSVD_Basis.npy', Mixture_foutname=today+'NC'+str(Nc)+'_NNDSVD_Mixture.npy')
