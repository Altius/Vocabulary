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

#from OONMF_afew_randos.py

list_of_realizations = [15,15,15,16,16,16,17,17,17,18,18,18]
list_of_seeds = [6,66,666]*4


for i,Nc in enumerate(list_of_realizations):

    a = OONMF.NMFobject(theNcomps=Nc)

    a.matrix_input_name( Basis_finname='92517/2017-09-28NC'+str(Nc)+'seed'+str(list_of_seeds[i])+'Basis.npy', Mixture_finname='92517/2017-09-28NC'+str(Nc)+'seed'+str(list_of_seeds[i])+'Mixture.npy')
    a.read_matrix_input()

    print(a.Basis.shape, 'basis shape')
    print(a.Mixture.shape, 'Mixture shape')

    a.Basis_Names = pd.read_table('colTitles.651samples_friendly.txt', header=None)[1].values

    basis_barsortorder = get_barsortorder(a.Basis)

    #print (len(barsortorder))

    a.make_stacked_bar_plot(a.BasisD, a.Basis.T, today+'Basis'+str(a.Ncomps)+'seed'+str(list_of_seeds[i])+'.png', barsortorder=basis_barsortorder, names=a.Basis_Names)

    normedBasis = (a.Basis.T / np.sum(a.Basis.T, axis=0)).T

    a.make_stacked_bar_plot(a.BasisD, normedBasis.T, today+'normedBasis'+str(a.Ncomps)+'seed'+str(list_of_seeds[i])+'.png', barsortorder=basis_barsortorder, names=a.Basis_Names)


    PRC = a.precision_recall_curve( A, names=a.Basis_Names, writefile=True, filename_addon='seed'+str(list_of_seeds[i]))

