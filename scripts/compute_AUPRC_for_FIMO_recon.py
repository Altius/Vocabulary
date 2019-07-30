# This script computes an area under precision recall curve for Motif Cluster Projection (formerly known as 2.5 D NMF) 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read in motif cluster occurence in masterlist
FIMO_ML = pd.read_csv('20619_FIMOtable_clusters.txt', sep='\t')

# load DHS mixture and motif projection mixture matrices 
Mixture = np.load('2018-06-08NC16_NNDSVD_Mixture.npy')
MotifMat = np.load('2019-05-08MotifClusterPA_X_ComponentsNC16NNDSVD.npy')

# calculate reconstruction
FIMO_RM = np.dot(MotifMat, Mixture).T

# since we're only doing presence absence, reduce the the clusters to a 1/0 
FIMO_ML_PA = FIMO_ML > 0.5 



# now use various decision boundaries to compute precision, recall for all motif clusters
                             
threshes_to_do = [0.001, 0.003, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99]

prec_table = []
rec_table = []
for thresh in threshes_to_do:
    print('doing ',thresh)
    PM05 = FIMO_RM > thresh
    
    TP_ar = np.multiply(PM05, FIMO_ML_PA.values)
    TN_ar = np.multiply(~PM05, ~FIMO_ML_PA.values.astype(bool))
    FP_ar = np.multiply(PM05, ~FIMO_ML_PA.values.astype(bool))
    FN_ar = np.multiply(~PM05, FIMO_ML_PA.values.astype(bool))
    
    TP_per_Motif = np.sum(TP_ar, axis=0) 
    FP_per_Motif = np.sum(FP_ar, axis=0) 
    TN_per_Motif = np.sum(TN_ar, axis=0) 
    FN_per_Motif = np.sum(FN_ar, axis=0) 
    
    prec_per_Motif = TP_per_Motif / (TP_per_Motif + FP_per_Motif + 1e-12)
    rec_per_Motif = TP_per_Motif / (TP_per_Motif + FN_per_Motif + 1e-12)
    
    prec_table.append(prec_per_Motif)
    rec_table.append(rec_per_Motif)
    
prec_table = np.array(prec_table)
rec_table = np.array(rec_table)

motif_PRAUC_ar = []
# integrate under precision recall curves
for i in range(prec_table.shape[1]):
    motif_PRAUC_ar.append(np.trapz([1] + list(rec_table[:,i]) + [0], [0] + list(prec_table[:,i]) +[1]))
    
# save results 

my_AUPRC_table = pd.DataFrame(np.array([FIMO_ML.columns.values, motif_PRAUC_ar]).T, columns=['Motif', 'AUPRC'])
my_AUPRC_table.to_csv('52019_motifAUPRC.csv', sep='\t')

