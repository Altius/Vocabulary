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

a = OONMF.NMFobject(theNcomps=16)

a.matrix_input_name(Basis_finname='2017-09-25NMF_Ncomps16Basis.npy', Mixture_finname='2017-09-25NMF_Ncomps16Mixture.npy')

a.read_matrix_input()

a.Basis_Names = pd.read_table('colTitles.651samples_friendly.txt', header=None)[1].values
a.normalize_matrices()

a.find_modules(A)
