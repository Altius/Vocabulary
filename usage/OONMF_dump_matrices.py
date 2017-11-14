import numpy as np
import pandas as pd
import sys 

import OONMF
from OONMFhelpers import *
today  = get_today()


a = OONMF.NMFobject(theNcomps=16)

a.matrix_input_name(Basis_finname='2017-09-25NMF_Ncomps16Basis.npy', Mixture_finname='2017-09-25NMF_Ncomps16Mixture.npy')
a.read_matrix_input()

a.Basis_Names = pd.read_table('colTitles.651samples_friendly.txt', header=None)[1].values


a.writeNMF_CSV(Basis_foutname='2017-09-25NMF_Ncomps16Basis.csv', Mixture_foutname='2017-09-25NMF_Ncomps16Mixture.csv')

a.normalize_matrices()
a.writeNMFnormed_CSV(Basis_foutname='2017-09-25NMF_Ncomps16Basis_Normed.csv', Mixture_foutname='2017-09-25NMF_Ncomps16Mixture_Normed.csv')

a.compute_reweighted_matrices()
a.writeNMFreweighted_CSV(Basis_foutname='2017-09-25NMF_Ncomps16Basis_Reweighted.csv', Mixture_foutname='2017-09-25NMF_Ncomps16Mixture_Reweighted.csv')
#may as well visualize these 

print('regular ',a.Basis.shape, a.Mixture.shape)
print('reweighted ',a.ReweightedBasis.shape, a.ReweightedMixture.shape)


basis_barsortorder = get_barsortorder(a.ReweightedBasis)

a.make_stacked_bar_plot(a.BasisD, a.ReweightedBasis.T, today+'ReweightedBasis'+str(a.Ncomps)+'.png', barsortorder=basis_barsortorder, names=a.Basis_Names)


a.normalize_reweighted_matrices()
a.writeNMFreweighted_normed_CSV(Basis_foutname='2017-09-25NMF_Ncomps16Basis_Reweighted_Normed.csv', Mixture_foutname='2017-09-25NMF_Ncomps16Mixture_Reweighted_Normed.csv')

a.make_stacked_bar_plot(a.BasisD, a.ReweightedNormedBasis.T, today+'ReweightedNormedBasis'+str(a.Ncomps)+'.png', barsortorder=basis_barsortorder, names=a.Basis_Names)

