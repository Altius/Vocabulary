'''
class: NMFobject

functions:
__init__
matrix_input_name
read_matrix_input
performNMF
build_reconstruction
normalize_matrices
normalize_reweighted_matrices
writeNMF_CSV
writeNMFnormed_CSV
writeNMFreweighted_CSV
writeNMFreweighted_normed_CSV
find_modules

'''

ClusterMode = True
import sys
#import datetime
#import time
import numpy as np
import pandas as pd
#from datetime import date
if (ClusterMode):
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.decomposition import NMF


from matplotlib.colors import LogNorm
from OONMFhelpers import *
today = get_today()


class NMFobject:
    def __init__(self, theNcomps=16):
        self.Basis = []
        self.Mixture = []
        self.Ncomps = theNcomps
        
        self.BasisD = 0
        self.MixtureD = 0
        
        self.Basis_Names = []
        self.Mixture_Names = []
        
        self.Reconstruction = []
        
        self.ReweightedBasis = []
        self.ReweightedMixture = []
        
        self.NormedBasis = []
        self.NormedMixture = []

        self.ReweightedNormedBasis = []
        self.ReweightedNormedMixture = []
    
    
    
    
    def matrix_input_name(self, Basis_finname='', Mixture_finname=''):
        if (len(Basis_finname) < 1 or len(Mixture_finname) < 1):
            #raise ValueError('syntax: read_matrix(Basis_finname, Mixture_finname)')
            print('syntax: read_matrix(Basis_finname, Mixture_finname)')
        self.Basis_finname = Basis_finname
        self.Mixture_finname = Mixture_finname

    def read_matrix_input(self):
        self.Basis = np.load(self.Basis_finname)
        self.BasisD = self.Basis.shape[0]
        if (len(self.Mixture_finname)>0):
            self.Mixture = np.load(self.Mixture_finname)
            self.MixtureD = self.Mixture.shape[1]
        
        
        
    def performNMF(self, data, randomseed=0):
        if(len(self.Basis) > 0):
            print('you are overwriting the Basis',self.Basis)
            cont = input('are you sure?')
            if (cont == 'n'):
                return
            
        model = NMF(n_components=self.Ncomps, init='random', random_state=randomseed)
        print('starting NMF at ', mytime(), flush=True)
        self.Basis = model.fit_transform(data) 
        print('done with NMF at ', mytime(), flush=True)
        self.Mixture = model.components_
        self.BasisD = self.Basis.shape[0]
        self.MixtureD = self.Mixture.shape[1]
        
    def build_reconstruction(self):
        self.Reconstruction = np.dot(self.Basis, self.Mixture)
                

    def normalize_matrices(self):
        self.NormedMixture =   self.Mixture / np.sum(self.Mixture, axis=0)
        self.NormedBasis =   (self.Basis.T / np.sum(self.Basis.T, axis=0)).T

    def normalize_reweighted_matrices(self):
        self.ReweightedNormedMixture =   self.ReweightedMixture / np.sum(self.ReweightedMixture, axis=0)
        self.ReweightedNormedBasis =   (self.ReweightedBasis.T / np.sum(self.ReweightedBasis.T, axis=0)).T


    def writeNMF(self, Basis_foutname, Mixture_foutname):
        np.save(Basis_foutname, self.Basis)
        #very confusing but it must be Mixture here for internal self-consisten\ cy. Can be Mixture.T for CSV files
        np.save(Mixture_foutname, self.Mixture)
        
        
        
    def writeNMF_CSV(self, Basis_foutname, Mixture_foutname):
        pd.DataFrame(self.Basis).to_csv(Basis_foutname)
        pd.DataFrame(self.Mixture.T).to_csv(Mixture_foutname)

    def writeNMFnormed_CSV(self, Basis_foutname, Mixture_foutname):
        pd.DataFrame(self.NormedBasis).to_csv(Basis_foutname)
        pd.DataFrame(self.NormedMixture.T).to_csv(Mixture_foutname)

    def writeNMFreweighted_CSV(self, Basis_foutname, Mixture_foutname):
        pd.DataFrame(self.ReweightedBasis).to_csv(Basis_foutname)
        pd.DataFrame(self.ReweightedMixture.T).to_csv(Mixture_foutname)

    def writeNMFreweighted_normed_CSV(self, Basis_foutname, Mixture_foutname):
        pd.DataFrame(self.ReweightedNormedBasis).to_csv(Basis_foutname)
        pd.DataFrame(self.ReweightedNormedMixture.T).to_csv(Mixture_foutname)


    def define_colorsA(self, mode='newSasha'):
        if (mode=='newSasha'):
            self.Comp_colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6','#6a3d9a','#ffff99','#b15928','#ffd700', '#AAAAAA', '#C52892', '#00bbbb']
        elif (mode=='Sasha'):
            self.Comp_colors = ['red', 'tan', 'lime','blue','m','k','c', 'coral', 'indigo','darkgreen','orange','grey','gold', 'lightskyblue', 'peru', 'olive']
        else:
            self.Comp_colors = ["#A6CEE3", "#438EC0", "#63A8A0", "#98D277", "#3BA432", "#B89B74", "#F16667", "#E62F27", "#F9A963", "#FE982C", "#ED8F47", "#C3AAD2", "#7D54A5","#B9A499", "#EAD27A" ,"#B15928"]
        if (self.Ncomps>16):
            np.random.seed(666)
            from matplotlib import colors as mcolors
            colornames = list(mcolors.CSS4_COLORS.keys())
            count = 16
            while (count < self.Ncomps):
                newcolor = colornames[np.random.randint(0,len(colornames))]
                trialcount = 0
                while ((newcolor in self.Comp_colors) and (trialcount < 100)):
                    newcolor = colornames[np.random.randint(0,len(colornames))]
                    trialcount+=1
                self.Comp_colors.append(newcolor)
                count+=1



    def make_stacked_bar_plot(self, Nrelevant, BarMatrix, bargraph_out, names = [], plotClusterMode=False, barsortorder=[], clusterTopLabels=[], colormode='newSasha'):
        if len(barsortorder)<1:
            barsortorder = np.arange(Nrelevant)
        if len(names) < 1:
            names = [str(i) for i in range(Nrelevant)]
            names = np.array(names)
        ttt = np.arange(Nrelevant)
        start = 0
        end = Nrelevant
        ground_pSample = ttt*0
        self.define_colorsA(mode=colormode)
        plt.clf()
        plt.figure(figsize=(150,40))
        plt.bar(ttt[start:end], BarMatrix[0,start:end][barsortorder], color=self.Comp_colors[0],
                 bottom=ground_pSample[start:end], alpha=1.0)
        ground_pSample = BarMatrix[0,start:end][barsortorder]
        for i in range(1,self.Ncomps):
            plt.bar(ttt[start:end],BarMatrix[i,start:end][barsortorder], bottom = ground_pSample, color=self.Comp_colors[i], alpha=1.0)
            ground_pSample = np.sum(BarMatrix[0: i+1,start:end], axis=0)[barsortorder]
        increase_axis_fontsize()
        plt.ylabel('sum of signal in matrix',fontsize=70)
        #plt.title('Full Sample',fontsize=70)
        samplenamesize = 11
        samplenamesize = (1/Nrelevant)**0.5 * 300
        #thebottom = 0.15
        thebottom = min([(1/Nrelevant)**0.3 * 1.2, 0.3])
        if(plotClusterMode):
            plt.xticks(ttt, ttt.astype(str), rotation='vertical', fontsize=samplenamesize)
            if len(clusterTopLabels) > 0:
                ax = plt.gca()
                ax2 = ax.twiny()
                ax2.set_xticks(ttt)
                ax2.set_xticklabels(clusterTopLabels.astype(str), rotation=90, fontsize=samplenamesize)
                #ax.xaxis.tick_top()
                #plt.xticks(ttt, clusterTopLabels.astype(str), rotation='vertical', fontsize=samplenamesize)
        else:
            plt.xticks(ttt, names[barsortorder], rotation='vertical', fontsize=samplenamesize)	
        plot_margin = 5
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=thebottom)
        plt.savefig(bargraph_out)
        plt.close()	
        
        
        
        
        
    def precision_recall_curve(self, data, names=[], writefile=False, filename_addon=''):
        #only works when objective matrix is known, and consists of entries 0/1. 
        if (len(self.Reconstruction)<1):
            self.build_reconstruction()
        print(data.shape, 'data')
        print(self.Reconstruction.shape, 'reconstruction')
        if (data.shape[0] != self.Reconstruction.shape[0] or data.shape[1] != self.Reconstruction.shape[1]):
            print('error! data and reconstruciton dont have matching shapes', data.shape, self.Reconstruction.shape )
            return
        if (np.max(data) > 1 or np.min(data) < 0):
            print('error, precision-recall curve only works for data between 0 and 1')
            return

        customthreshes = [0.05, 0.1, 0.15, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
        recall_ar = []
        precision_ar = []

        if len(names) < 1:
            print('filling in names')
            names = np.arange(self.BasisD).astype(str)
        prec_recall_table = []
        for customthresh in customthreshes:
            F1_ar = []
            if(writefile):
                f = open(today+'SampleCustomStats_thresh'+str(customthresh)+'_Ncomp'+str(self.Ncomps)+filename_addon+'.txt', 'w')
            count = 0

            for sample in range(self.BasisD):
                DHSar_cut = data[sample]>customthresh
                predDHSar_cut = self.Reconstruction[sample] > customthresh


                TP = len(self.Reconstruction[sample][DHSar_cut * predDHSar_cut]) 
                FP = len(self.Reconstruction[sample][predDHSar_cut * np.invert( DHSar_cut)])
                TN = len(self.Reconstruction[sample][np.invert(predDHSar_cut) * np.invert(DHSar_cut)])
                FN = len(self.Reconstruction[sample][np.invert(predDHSar_cut) * (DHSar_cut)])

                if ((TP + FN) > 0 ):
                    recall = TP / (TP + FN)
                else: 
                    recall = 0
                if ((TP + FP) > 0 ):
                    precision = TP / (TP + FP)
                else:
                    precision=0
                accuracy = (TP + TN) /(len(self.Reconstruction[sample]))
                if (precision + recall) == 0:
                    F1 = 0
                else:
                    F1 = 2*(precision*recall)/(precision+recall)
                F1_ar.append(F1)
                if (writefile):
                    print(sample, names[count], TP, FP, TN, FN, recall, precision, accuracy, F1, file=f)
                prec_recall_table.append([customthresh, sample, names[count], TP, FP, TN, FN, recall, precision, accuracy, F1])
                count +=1
            print('Ncomps ',self.Ncomps, 'thresh ', customthresh, ' mean F1 score ',np.mean(np.array(F1_ar)))
        return pd.DataFrame(prec_recall_table, columns=['threshold', 'sample_number', 'sample_name', 'TP', 'FP', 'TN', 'FN', 'recall', 'precision', 'accuracy', 'F1'])
        
        
        
        
    def compute_reweighted_matrices(self):
    
        bigAllDHSSum_ar = []
        bigAllSampleSum_ar = []
        
        #formerly known as sansvar. The point of this is to attribute "new meaning" to components by reweighting the Basis by the Mixture and vice versa. Doing this allows for an interpretation where each component-vector in the Basis reflects #DHSs found in that sample and vice versa. 
        
        for i in range(self.Ncomps):
            bongo = np.copy(self.Basis)
            for k in range(self.Ncomps):
                if (k!=i):
                    bongo[:,k]*=0
            sansvar = np.dot(bongo, self.Mixture)
            bigAllDHSSum_ar.append(np.sum(sansvar[:,0:], axis=1))
            bigAllSampleSum_ar.append(np.sum(sansvar[:,0:], axis=0))
            del(sansvar)
        
        self.ReweightedBasis = np.array(bigAllDHSSum_ar).T
        self.ReweightedMixture = np.array(bigAllSampleSum_ar)
            


    def find_modules(self, data, ClustMult=4, chosenthresh = 0.35, cosdist_sample_thresh=0.3, cosdist_DHS_thresh=0.3 ):
    
        #initially only implementing SampleNormed basis
        #but with a fix to make the DHSappearCut  more self-consistent
        
        from sklearn.cluster import KMeans
        import scipy.spatial.distance as spdist

        kmeans_normed_sample_orig = KMeans(n_clusters=self.Ncomps*ClustMult, random_state=0).fit(self.NormedBasis)
        ClusterMeans = kmeans_normed_sample_orig.cluster_centers_
        #this doesn't really work... i am using normed versionf of vectors and then throwing them into DHS matrix for reconstruction
        ClusterPredictions = kmeans_normed_sample_orig.predict(relevant_matrix)
        
        Clusters = []
        # find the "cluster center" in the actual DHS decomposed space
        ClusterMeansReal = []

        for i in range(self.Ncomps*ClustMult):
            current_cut = np.argwhere(ClusterPredictions==i).T[0]
            #print (BasisMat[current_cut].shape)
            themean = np.mean(self.Basis[current_cut], axis=0)
            #print(themean.shape)
            print(i, len(ClusterPredictions[current_cut]), names[current_cut])
            ClusterMeansReal.append(themean)
            Clusters.append([i, self.Basis_Names[current_cut], current_cut ])        
        
        #using ClusterMeans  (from Kmeans)

        #threshold for what counts as a positive prediction in reconstruction matrix
        #chosenthresh = 0.35
        #threshold for cosine similarity distance for samples to be included in this module. 0.1 gives about the same samples as in the original cluster for 80% component dominance. 
        #cosdist_sample_thresh = 0.3

        #threshold for cosine similarity distance for DHS to be included in this module
        #cosdist_DHS_thresh = 0.3
        for i, cluster in enumerate(Clusters):

            print ('doing cluster ',i)
            # mean sample-wise component vector for this cluster
            mean_of_cluster = ClusterMeans[i]
            real_mean_of_cluster = ClusterMeansReal[i]
            #create a single 1 x NDHS prediction for the cluster to
            recon_cluster_DHS = np.dot(mean_of_cluster, self.Mixture)
            cluster_length = len(cluster[1])
            print('cluster length ',cluster_length)
            print('name of samples', cluster[1])

            #now you tile the vector to a size equivalent to N_samples_in_cluster x NDHS
            kar = np.tile(recon_cluster_DHS, cluster_length).reshape(cluster_length, self.MixtureD)

            #***** NOTE: for some reason I decided to use original version of NMF matrices for cosine 
            #***** Similarity scores. Self-consistent? Not sure! 

            #calculate cosine similarity of cluster mean to all sample-component vectors 
            sample_cosdists = spdist.cdist(self.NormedBasis, np.reshape(mean_of_cluster, (1,self.Ncomps)), 'cosine').flatten()    
            #calculate cosine similarity of cluster mean to all DHS-component vectors 
            DHScosdists = spdist.cdist(self.NormedMixture, np.reshape(mean_of_cluster, (1,self.Ncomps)), 'cosine').flatten()
            #now we sort these sample cosine distances
            closers = np.argsort(sample_cosdists)

            #using the thresholds defined outside the loop, we find the samples that meet the closeness
            #criteria 
            Num_samples = len(sample_cosdists[sample_cosdists<cosdist_sample_thresh])
            print(' Num samples', Num_samples)
            print('samples meeting threshold ',names[closers[0:Num_samples]])

            #First we use our decision threshold to figure out if the model predicts a positive
            #for a given position. Using the sample cluster mean component vector
            DHSappear_cut = np.dot(real_mean_of_cluster, self.Mixture)>chosenthresh
            #now sort the coside distances of DHS's that will appear from previous step    
            DHScloseness = np.argsort(DHScosdists[DHSappear_cut])
            Num_DHSs = len(DHScosdists[DHSappear_cut][DHScosdists[DHSappear_cut] < cosdist_DHS_thresh])
            print(' Num_DHSs', Num_DHSs)

            #now we use the masks to create a miniature version of the reconstruction matrix for the given DHS only
            #strange use of double transpose 
            sorted_Module = np.dot(self.Basis[closers[0:Num_samples]], self.Mixture.T[DHSappear_cut][DHScloseness][0:Num_DHSs].T)
            #no idea how i pulled out this magic
            #this gives the real data corresponding to the module 
            corresponding_data = data[closers[0:Num_samples]][:,np.argwhere(DHSappear_cut == True).T[0][DHScloseness][0:Num_DHSs]]
            #flatten both modules to count accurate predictions. 
            flatsorted_Module = sorted_Module.flatten()
            flatcorresponding_data = corresponding_data.flatten()
            #calculate "density" i.e. number of entries of DHS hits within the module
            if (len(flatsorted_Module) > 0):
                density_in_recon = len(flatsorted_Module[flatsorted_Module>chosenthresh])/len(flatsorted_Module)
            else:
                density_in_recon = 0
            if(len(flatcorresponding_data) > 0):
                density_in_data = len(flatcorresponding_data[flatcorresponding_data>chosenthresh])/len(flatcorresponding_data)
            else:
                density_in_data = 0
            #calculate true positives, false positives, False negatives, precision ,recall
            TP = len(flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
            FP = len(flatcorresponding_data[ np.invert(flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
            FN = len( flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * np.invert(flatsorted_Module>chosenthresh)])
            if ((TP + FN ) > 0):
                recall = TP / (TP + FN)
            else: 
                recall = 0
            if ((TP +FP) > 0):
                precision = TP / (TP + FP)
            else:
                precision = 0
            print('density in reconstruction,  density in data,  precision,  recall')
            print(density_in_recon, density_in_data, precision, recall  )
            # don't remmeber how i figured this out either: I was possessed 
            #saves module as a list of samples in the header followed by a list of eahc DHS position
            #and whether it is or isn't present in the module.
            firstcut = np.copy(DHSappear_cut)
            firstcut[DHSappear_cut] *= (DHScosdists[DHSappear_cut] < cosdist_DHS_thresh)
            list1 = [Num_samples, Num_DHSs, chosenthresh, cosdist_sample_thresh,  cosdist_DHS_thresh, density_in_recon, density_in_data, precision, recall]
            theheader =  ' '.join(str(e) for e in list1)
            namestr = ' '.join(names[closers[0:Num_samples]])
            theheader +=' '+namestr
            foutname = today+'Cluster_SampleNormed_Module'+str(i)+'.txt'
            np.savetxt(foutname, firstcut, fmt='%i', header=theheader)

        foutname2= today+'SampleNormed_ModuleMeans_ClusterMeansReal.txt'
        np.savetxt(foutname2, np.array(ClusterMeansReal))
        foutname2= today+'SampleNormed_ModuleMeans_ClusterMeans.txt'
        np.savetxt(foutname2, np.array(ClusterMeans))
        
