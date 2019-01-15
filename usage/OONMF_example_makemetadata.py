# use OONMFmetadata library to make our signature category statistical association plots

import numpy as np
import pandas as pd
import sys

import OONMF
import OONMFhelpers as OH
import OONMFmetadata as OMD

today = OH.get_today()

myNC=16
usingCut = False
thecut = []

if len(sys.argv) > 1:
    finname_base = sys.argv[1]
    myNC=int(sys.argv[2])

if len(sys.argv)> 3:
    print('adding cut')
    usingCut = True
    cut_finname = sys.argv[3]
    thecut = np.load(cut_finname)
    
    
myNMFcomps=myNC

finname = 'ENCODE733MetaData.csv'
MetaData  = pd.read_csv(finname, sep='\t')
print (MetaData.head())

a = OONMF.NMFobject(theNcomps=myNMFcomps)
sampnamePD = pd.read_table('sampnams_733.txt', header=None, names=['LN', 'DS', 'type'])
sampnamePD['full_name'] = sampnamePD.LN + '-' + sampnamePD.DS + '-' + sampnamePD.type

thenames = sampnamePD.full_name.values
fullBasis_Names = thenames
if (usingCut):
    print('adjusting names')
    thenames = thenames[thecut]
a.Basis_Names = thenames
rosetta = OMD.get_rosetta(MetaData, fullBasis_Names)


print(a.Basis_Names[0:5])


foutname_base = finname_base+today+'_'  

a.matrix_input_name(Basis_finname=finname_base+'Basis.npy', Mixture_finname=finname_base+'Mixture.npy')

a.read_matrix_input()


if (usingCut):
    print('a' ,MetaData.shape, len(rosetta), len(thecut))
    newMetaData = MetaData.loc[rosetta][thecut]
    print('b', len(newMetaData))
    newrosetta = np.arange(len(newMetaData)).astype(int)

else:
    newMetaData = MetaData
    newrosetta = rosetta

# all categories that have proven interesting, and then some
categories = ['system', 'organ','Biological_state','subsystem','Sex','germ layer', 'Growth stage', 'class', 'lib_kit_method', 'Sample_group', 'Candidate removal -- insert sizes & dupe rate', 'Candidate removal -- insert sizes', 'Donor_ID', 'Ethnicity', 'lib_cleanup', 'sample_dhs_protease_inhibitor']

for my_cat_str in categories:
    MetaData[my_cat_str].fillna('None', inplace=True)
    (SystemChart,SystemTypes) = OMD.make_significance_plot(newMetaData, a.Basis, my_cat_str, newrosetta,thecmap='inferno', NMFCOMPS=myNMFcomps,filename_addon=foutname_base, write_mode=True)
    
    # at the moment, homogeneity is not very useful / defunct 
    #make_significance_plot_homogeneity(MetaData, homogeneity, my_cat_str, rosetta, thecmap='inferno', NMFCOMPS=myNMFcomps)
