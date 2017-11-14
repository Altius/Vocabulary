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
for i in range (2, 25):

	a = OONMF.NMFobject(theNcomps=i)

	a.matrix_input_name(Basis_finname='2017-09-25NMF_Ncomps'+str(i)+'Basis.npy', Mixture_finname='2017-09-25NMF_Ncomps'+str(i)+'Mixture.npy')
    
	a.read_matrix_input()

	a.Basis_Names = pd.read_table('colTitles.651samples_friendly.txt', header=None)[1].values

	basis_barsortorder = get_barsortorder(a.Basis)

	#print (len(barsortorder))

	a.make_stacked_bar_plot(a.BasisD, a.Basis.T, today+'Basis'+str(a.Ncomps)+'.png', barsortorder=basis_barsortorder, names=a.Basis_Names)

	normedBasis = (a.Basis.T / np.sum(a.Basis.T, axis=0)).T

	a.make_stacked_bar_plot(a.BasisD, normedBasis.T, today+'normedBasis'+str(a.Ncomps)+'.png', barsortorder=basis_barsortorder, names=a.Basis_Names)


	PRC = a.precision_recall_curve( A, names=Basis_names, writefile=True)

