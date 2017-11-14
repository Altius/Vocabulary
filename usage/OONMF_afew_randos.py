import numpy as np
import pandas as pd
import OONMF
from OONMFhelpers import *

today  = get_today()

def acquire_data():
    # fill in for whatever application
	pass
	print('starting read at ', mytime(), flush=True)
	A = pd.read_table('masterListPresenceAbsenceMatrix.651samples.FDR0.0010.hg38.bed', header=None)
	print('finished read at ', mytime(), flush=True)
	return A.drop([A.columns[0], A.columns[1], A.columns[2]], axis=1).values.T

A = acquire_data()
list_of_realizations = [15,15,15,16,16,16,17,17,17,18,18,18]
list_of_seeds = [6,66,666]*4
for i,Nc in enumerate(list_of_realizations):
	print(i)
	a = OONMF.NMFobject(theNcomps=Nc)
	a.performNMF(data=A, randomseed=list_of_seeds[i])	
	a.writeNMF(Basis_foutname= today+'NC'+str(Nc)+'seed'+str(list_of_seeds[i])+'Basis.npy', Mixture_foutname=today+'NC'+str(Nc)+'seed'+str(list_of_seeds[i])+'Mixture.npy')

