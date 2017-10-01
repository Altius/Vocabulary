import numpy as np
import pandas as pd

import OONMF
from OONMFhelpers import *
from OONMFmetadata import *

myNMFcomps=17

finname = 'ENCODE 651-sample list WM20170918 - ENCODE3 metadata.csv'
MetaData  = pd.read_csv(finname)
MetaData = MetaData[0:827]


a = OONMF.NMFobject(theNcomps=myNMFcomps)
a.Basis_Names = pd.read_table('colTitles.651samples_friendly.txt', header=None)[1].values

a.matrix_input_name(Basis_finname='92517/2017-09-25NMF_Ncomps'+str(myNMFcomps)+'Basis.npy', Mixture_finname='92517/2017-09-25NMF_Ncomps'+str(myNMFcomps)+'Mixture.npy')

a.read_matrix_input()

homogeneity = []
for sample in a.Basis:
    hom = np.max(sample)/np.sum(sample)
    homogeneity.append(hom)
homogeneity = np.array(homogeneity)

rosetta = get_rosetta(MetaData, a.Basis_Names)

categories = ['system', 'organ','Biological_state','subsystem','Sex','germ layer', 'Growth stage']

for my_cat_str in categories:
    MetaData[my_cat_str].fillna('None', inplace=True)
    (SystemChart,SystemTypes) = make_significance_plot(MetaData, a.Basis, my_cat_str, rosetta,thecmap='inferno', NMFCOMPS=myNMFcomps)
    make_significance_plot_homogeneity(MetaData, homogeneity, my_cat_str, rosetta, thecmap='inferno', NMFCOMPS=myNMFcomps)
