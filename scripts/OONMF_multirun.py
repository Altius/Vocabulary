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

#number of components (k) - list of values to try
list_of_realizations = [5,6,7,8]

#random seed (not very important for NNDSVD)
list_of_seeds = [20]*len(list_of_realizations)

#loop of runs
for i,Nc in enumerate(list_of_realizations):
	print(i)
	a = OONMF.NMFobject(theNcomps=Nc)
	a.performNMF(data=A, randomseed=list_of_seeds[i], theinit='nndsvd')	
	a.writeNMF(Basis_foutname= today+'NC'+str(Nc)+'_NNDSVD_Basis.npy', Mixture_foutname=today+'NC'+str(Nc)+'_NNDSVD_Mixture.npy')

