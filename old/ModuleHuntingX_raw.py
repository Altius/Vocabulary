# # first re-capitulate how to find modules with clusters of samples 

import sys
import datetime
import time



from sklearn.decomposition import NMF
import numpy as np
import pandas as pd
from datetime import date
today = str(date.today())

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

BelugaMode = False
#saveDHSmat = True


# several parameters
axis_fontsize = 30
myrandom = 0

def mytime():
	ts = time.time()
	sttime = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H:%M:%S')
	return sttime
print('starting at ', mytime(), flush=True)

def increase_axis_fontsize():
	ax = plt.gca()
	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(axis_fontsize)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(axis_fontsize)
		label.set_family('serif')
		
def make_bar_plot(Nrelevant, BarMatrix, bargraph_out, plotClusterMode=False, barsortorder=[], clusterTopLabels=[], myNMFcomps=16):
	if len(barsortorder)<1:
		barsortorder = np.arange(Nrelevant)
	ttt = np.arange(Nrelevant)
	start = 0
	end = Nrelevant
	ground_pSample = ttt*0
	plt.clf()
	plt.figure(figsize=(150,40))
	plt.bar(ttt[start:end], BarMatrix[0,start:end][barsortorder], color='r',
             bottom=ground_pSample[start:end], alpha=0.75)
	ground_pSample = BarMatrix[0,start:end][barsortorder]
	for i in range(1,myNMFcomps):
		plt.bar(ttt[start:end],BarMatrix[i,start:end][barsortorder], bottom = ground_pSample, color=Comp_colors[i], alpha=0.75)
		ground_pSample = np.sum(BarMatrix[0: i+1,start:end], axis=0)[barsortorder]
	increase_axis_fontsize()
	plt.ylabel('sum of signal in matrix',fontsize=70)
	#plt.title('Full Sample',fontsize=70)
	samplenamesize = 27
	thebottom = 0.25
	if (BelugaMode):
		samplenamesize = 11
		thebottom = 0.15
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

#Nsamples
if (BelugaMode):
	Nsamples = 827
	NMFcomps = 0 #will be assigned later
	B = pd.read_table('colTitles_friendly.txt', header=None)
	names = B[1].values
else:
	Nsamples = 176
	NMFcomps = 16 #has to be set here since NMF will be done on the spot
	B = pd.read_table('rowTitles.txt', header=None)
	names = B[0].values


if (BelugaMode):
	if (len(sys.argv) < 3):
		print ('syntax: blah.py pickle1(daslog) pickle2(daslog2)')
		sys.exit()
	finname1 = str(sys.argv[1])
	finname2 = str(sys.argv[2])
	print('starting read at ', mytime(), flush=True)
	A = pd.read_table('masterListPresenceAbsence_827samples_FDR0.0010_hg38.bed')
	print('finished read at ', mytime(), flush=True)
	print(A.head())
	print(A.shape)
	a = A.drop([A.columns[0], A.columns[1], A.columns[2]], axis=1).values
	a = a.T
	print('finished drop at ', mytime(), flush=True)
	del(A)
else:
	IslandSize = 400
	finname = 'mtxlist.txt'
	f = open(finname, 'r')
	dars = f.readlines()
	featurelist = []
	for line in dars:
		g = np.loadtxt(line.strip())
		#daslog = np.log10(g+1)
		daslog = (g > 0).astype(int) 
		if (daslog.shape[1] != IslandSize):
			daslog = daslog[:,:IslandSize]
		featurelist.append(daslog)
	featurelist = np.array(featurelist)
	a = np.copy(featurelist[0])
	for i in range(499):
		a= np.append(a, featurelist[i+1], axis=1)
	#print (a.shape)


if (BelugaMode):
	print('starting read NMF results at ', mytime(), flush=True)
	BasisMat = np.load(finname1)
	MixtureMat = np.load(finname2)
else:
	model = NMF(n_components=NMFcomps, init='random', random_state=0)
	BasisMat = model.fit_transform(a) 
	MixtureMat = model.components_
	print(model.reconstruction_err_)

reconstruction = np.dot(BasisMat, MixtureMat)
NMFcomps = BasisMat.shape[1]
NDHS = MixtureMat.shape[1]

Comp_colors = ['red', 'tan', 'lime','blue','m','k','c', 'coral', 'indigo','darkgreen','orange','grey','gold', 'lightskyblue', 'peru', 'olive']
if (NMFcomps>16):
	from matplotlib import colors as mcolors
	colornames = list(mcolors.CSS4_COLORS.keys())
	count = 16
	while (count < NMFcomps):
		newcolor = colornames[np.random.randint(0,len(colornames))]
		trialcount = 0
		while ((newcolor in Comp_colors) and (trialcount < 100)):
			newcolor = colornames[np.random.randint(0,len(colornames))]
			trialcount+=1
		Comp_colors.append(newcolor)
		count+=1
		
#for now lets just do unweighted version


# In[2]:


#calculate normed version of Sample-wise component vectors
rawnormed_Basisums = np.sum(BasisMat, axis=1)
rawnormed_VecRats = np.copy(BasisMat)
for i, rat in  enumerate(rawnormed_VecRats):
	if (rawnormed_Basisums[i]>0):
		rat/=rawnormed_Basisums[i]

#calculate normed version of DHS-wise component vectors
rawnormed_Mixturesums = np.sum(MixtureMat.T, axis=1)
rawnormed_MixtureRats = np.copy(MixtureMat.T)
rawnormed_MixtureRats = rawnormed_MixtureRats[rawnormed_Mixturesums>0]
rawnormed_Mixturesums = rawnormed_Mixturesums[rawnormed_Mixturesums>0]
for i, rat in  enumerate(rawnormed_MixtureRats):
    rat/=rawnormed_Mixturesums[i]

from sklearn.cluster import KMeans

kmeans_normed_sample_orig = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(rawnormed_VecRats)
ClusterMeans = kmeans_normed_sample_orig.cluster_centers_
#this doesn't really work... i am using normed versionf of vectors and then throwing them into DHS matrix for reconstruction

ClusterPredictions = kmeans_normed_sample_orig.predict(rawnormed_VecRats)



# In[39]:

#which things in which cluster

Clusters = []
# find the "cluster center" in the actual DHS decomposed space
ClusterMeansReal = []

for i in range(NMFcomps*2):
    current_cut = np.argwhere(ClusterPredictions==i).T[0]
    #print (BasisMat[current_cut].shape)
    themean = np.mean(BasisMat[current_cut], axis=0)
    #print(themean.shape)
    print(i, len(ClusterPredictions[current_cut]), names[current_cut])
    ClusterMeansReal.append(themean)
    Clusters.append([i, names[current_cut], current_cut ])



# In[13]:

import scipy.spatial.distance as spdist


# In[19]:

# In[40]:

#using ClusterMeans Real (actual values from newdaslog)

#threshold for what counts as a positive prediction in reconstruction matrix
chosenthresh = 0.5
#threshold for cosine similarity distance for samples to be included in this module. 0.1 gives about the same samples as in the original cluster for 80% component dominance. 
cosdist_sample_thresh = 0.3

#threshold for cosine similarity distance for DHS to be included in this module
cosdist_DHS_thresh = 0.3
for i, cluster in enumerate(Clusters):

	print ('doing cluster ',i)
	# mean sample-wise component vector for this cluster
	mean_of_cluster = ClusterMeansReal[i]
	#create a single 1 x NDHS prediction for the cluster to
	recon_cluster_DHS = np.dot(mean_of_cluster, MixtureMat)
	cluster_length = len(cluster[1])
	print('cluster length ',cluster_length)
	print('name of samples', cluster[1])
    
	#now you tile the vector to a size equivalent to N_samples_in_cluster x NDHS
	kar = np.tile(recon_cluster_DHS, cluster_length).reshape(cluster_length, NDHS)

	#***** NOTE: for some reason I decided to use original version of NMF matrices for cosine 
	#***** Similarity scores. Self-consistent? Not sure! 

	#calculate cosine similarity of cluster mean to all sample-component vectors 
	sample_cosdists = spdist.cdist(BasisMat, np.reshape(mean_of_cluster, (1,NMFcomps)), 'cosine').flatten()    
	#calculate cosine similarity of cluster mean to all DHS-component vectors 
	DHScosdists = spdist.cdist(MixtureMat.T, np.reshape(mean_of_cluster, (1,NMFcomps)), 'cosine').flatten()
	#now we sort these sample cosine distances
	closers = np.argsort(sample_cosdists)

	#using the thresholds defined outside the loop, we find the samples that meet the closeness
	#criteria 
	Num_samples = len(sample_cosdists[sample_cosdists<cosdist_sample_thresh])
	print(' Num samples', Num_samples)
	print('samples meeting threshold ',names[closers[0:Num_samples]])

	#First we use our decision threshold to figure out if the model predicts a positive
	#for a given position. Using the sample cluster mean component vector
	DHSappear_cut = np.dot(mean_of_cluster, MixtureMat)>chosenthresh
	#now sort the coside distances of DHS's that will appear from previous step    
	DHScloseness = np.argsort(DHScosdists[DHSappear_cut])
	Num_DHSs = len(DHScosdists[DHSappear_cut][DHScosdists[DHSappear_cut] < cosdist_DHS_thresh])
	print(' Num_DHSs', Num_DHSs)

	#now we use the masks to create a miniature version of the reconstruction matrix for the given DHS only
	#strange use of double transpose 
	sorted_Module = np.dot(BasisMat[closers[0:Num_samples]], MixtureMat.T[DHSappear_cut][DHScloseness][0:Num_DHSs].T)
	#no idea how i pulled out this magic
	#this gives the real data corresponding to the module 
	corresponding_data = a[closers[0:Num_samples]][:,np.argwhere(DHSappear_cut == True).T[0][DHScloseness][0:Num_DHSs]]
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
	TP =len(flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
	FP =len(flatcorresponding_data[ np.invert(flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
	FN =len(flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * np.invert(flatsorted_Module>chosenthresh)])
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
	foutname = today+'Cluster_SampleRaw_Module'+str(i)+'.txt'
	np.savetxt(foutname, firstcut, fmt='%i', header=theheader)


# # now find modules starting from DHS clusters



# In[21]:

kmeans_normed_DHS_orig = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(rawnormed_MixtureRats)


# In[22]:

ClusterMeansDHS = kmeans_normed_DHS_orig.cluster_centers_
ClusterPredictionsDHS = kmeans_normed_DHS_orig.predict(rawnormed_MixtureRats)



ClustersDHS = []
ClusterMeanRealDHS = []
for i in range(NMFcomps*2):
    current_cut = np.argwhere(ClusterPredictionsDHS==i).T[0]
    themean = np.mean(MixtureMat.T[current_cut], axis=0)
    ClusterMeanRealDHS.append(themean)
    print(len(current_cut), themean.shape)
    print(i, len(ClusterPredictionsDHS[current_cut]), current_cut)
    #to mimic sample names, I'm just using the second element to "name" DHS's
    ClustersDHS.append([i, current_cut.astype(str), current_cut])




# # trying again via loop (using real DHS centroids) 

# In[133]:

chosenthresh = 0.15
#threshold for cosine similarity distance for samples to be included in this module. 0.1 gives about the same samples as in the original cluster for 80% component dominance. 
cosdist_sample_thresh = 0.75

#threshold for cosine similarity distance for DHS to be included in this module
cosdist_DHS_thresh = 0.3
for i, cluster in enumerate(ClustersDHS):

	print ('doing cluster ',i)
	# mean sample-wise component vector for this cluster
	mean_of_cluster = ClusterMeanRealDHS[i]
	#create a single 1 x NDHS prediction for the cluster to
	recon_cluster_Sam = np.dot(BasisMat, mean_of_cluster)
	#recon_cluster_DHS = np.dot(mean_of_cluster, MixtureMat)

	cluster_length = len(cluster[1])
	print('cluster length ',cluster_length)
	print('name of samples', cluster[1])

	#now you tile the vector to a size equivalent to N_DHS_in_cluster x NSamples
	kar = np.tile(recon_cluster_Sam, cluster_length).reshape(cluster_length, Nsamples)
	print(kar.shape)
	kar = kar.T

	#print(kar[0]-kar[1])    
	#***** NOTE: for some reason I decided to use original version of NMF matrices for cosine 
	#***** Similarity scores. Self-consistent? Not sure! 

	#calculate cosine similarity of cluster mean to all sample-component vectors 
	sample_cosdists = spdist.cdist(BasisMat, np.reshape(mean_of_cluster, (1,NMFcomps)), 'cosine').flatten()    
	#calculate cosine similarity of cluster mean to all DHS-component vectors 
	DHScosdists = spdist.cdist(MixtureMat.T, np.reshape(mean_of_cluster, (1,NMFcomps)), 'cosine').flatten()
	#now we sort these sample cosine distances

	closers = np.argsort(DHScosdists)

	#using the thresholds defined outside the loop, we find the samples that meet the closeness
	#criteria 
	Num_DHSs = len(DHScosdists[DHScosdists<cosdist_sample_thresh])
	print(' Num DHSs', Num_DHSs)
	print('DHSs meeting threshold ',closers[0:Num_DHSs])
	#First we use our decision threshold to figure out if the model predicts a positive
	#for a given position. Using the sample cluster mean component vector
	Sample_appear_cut = np.dot(mean_of_cluster, BasisMat.T)>chosenthresh
	#now sort the coside distances of DHS's that will appear from previous step    

	Sample_closeness = np.argsort(sample_cosdists[Sample_appear_cut])
	print('close samples ',sample_cosdists[Sample_appear_cut])
	Num_samples = len(sample_cosdists[Sample_appear_cut][sample_cosdists[Sample_appear_cut] < cosdist_sample_thresh])
	print(' Num_samples', Num_samples)

	sorted_Module = np.dot(BasisMat[Sample_appear_cut][Sample_closeness][0:Num_samples], MixtureMat.T[closers][0:Num_DHSs].T)
	corresponding_data = a[Sample_appear_cut][Sample_closeness][0:Num_samples][:,np.arange(NDHS)[closers][0:Num_DHSs]]
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
	TP =len(flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
	FP =len(flatcorresponding_data[ np.invert(flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
	FN =len(flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * np.invert(flatsorted_Module>chosenthresh)])
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
	firstcut = np.zeros(NDHS)
	firstcut[closers][0:Num_DHSs] += 1
	list1 = [Num_samples, Num_DHSs, chosenthresh, cosdist_sample_thresh,  cosdist_DHS_thresh, density_in_recon, density_in_data, precision, recall]
	theheader =  ' '.join(str(e) for e in list1)
	namestr = ' '.join(names[Sample_appear_cut][Sample_closeness][0:Num_samples])
	theheader +=' '+namestr
	print(namestr)
	foutname = today+'Cluster_DHSRaw_Module'+str(i)+'.txt'
	np.savetxt(foutname, firstcut, fmt='%i', header=theheader)


# # now even though it is in consistent, lets use the normalized DHS vectors for cluster means

# In[134]:

chosenthresh = 0.3
#threshold for cosine similarity distance for samples to be included in this module. 0.1 gives about the same samples as in the original cluster for 80% component dominance. 
cosdist_sample_thresh = 0.5

#threshold for cosine similarity distance for DHS to be included in this module
cosdist_DHS_thresh = 0.3
for i, cluster in enumerate(ClustersDHS):

	print ('doing cluster ',i)
	# mean sample-wise component vector for this cluster
	mean_of_cluster = ClusterMeansDHS[i]
	#create a single 1 x NDHS prediction for the cluster to
	recon_cluster_Sam = np.dot(BasisMat, mean_of_cluster)
	#recon_cluster_DHS = np.dot(mean_of_cluster, MixtureMat)

	cluster_length = len(cluster[1])
	print('cluster length ',cluster_length)
	print('name of samples', cluster[1])

	#now you tile the vector to a size equivalent to N_DHS_in_cluster x NSamples
	kar = np.tile(recon_cluster_Sam, cluster_length).reshape(cluster_length, Nsamples)
	print(kar.shape)
	kar = kar.T

	#print(kar[0]-kar[1])    
	#***** NOTE: for some reason I decided to use original version of NMF matrices for cosine 
	#***** Similarity scores. Self-consistent? Not sure! 

	#calculate cosine similarity of cluster mean to all sample-component vectors 
	sample_cosdists = spdist.cdist(BasisMat, np.reshape(mean_of_cluster, (1,NMFcomps)), 'cosine').flatten()    
	#calculate cosine similarity of cluster mean to all DHS-component vectors 
	DHScosdists = spdist.cdist(MixtureMat.T, np.reshape(mean_of_cluster, (1,NMFcomps)), 'cosine').flatten()
	#now we sort these sample cosine distances

	closers = np.argsort(DHScosdists)

	#using the thresholds defined outside the loop, we find the samples that meet the closeness
	#criteria 
	Num_DHSs = len(DHScosdists[DHScosdists<cosdist_sample_thresh])
	print(' Num DHSs', Num_DHSs)
	print('DHSs meeting threshold ',closers[0:Num_DHSs])
	#First we use our decision threshold to figure out if the model predicts a positive
	#for a given position. Using the sample cluster mean component vector
	Sample_appear_cut = np.dot(mean_of_cluster, BasisMat.T)>chosenthresh
	#now sort the coside distances of DHS's that will appear from previous step    

	Sample_closeness = np.argsort(sample_cosdists[Sample_appear_cut])
	print('close samples ',sample_cosdists[Sample_appear_cut])
	Num_samples = len(sample_cosdists[Sample_appear_cut][sample_cosdists[Sample_appear_cut] < cosdist_sample_thresh])
	print(' Num_samples', Num_samples)

	sorted_Module = np.dot(BasisMat[Sample_appear_cut][Sample_closeness][0:Num_samples], MixtureMat.T[closers][0:Num_DHSs].T)
	corresponding_data = a[Sample_appear_cut][Sample_closeness][0:Num_samples][:,np.arange(NDHS)[closers][0:Num_DHSs]]
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
	TP =len(flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
	FP =len(flatcorresponding_data[ np.invert(flatcorresponding_data>chosenthresh) * (flatsorted_Module>chosenthresh)])
	FN =len(flatcorresponding_data[ (flatcorresponding_data>chosenthresh) * np.invert(flatsorted_Module>chosenthresh)])
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
	firstcut = np.zeros(NDHS)
	firstcut[closers][0:Num_DHSs] += 1
	list1 = [Num_samples, Num_DHSs, chosenthresh, cosdist_sample_thresh,  cosdist_DHS_thresh, density_in_recon, density_in_data, precision, recall]
	theheader =  ' '.join(str(e) for e in list1)
	namestr = ' '.join(names[Sample_appear_cut][Sample_closeness][0:Num_samples])
	theheader +=' '+namestr
	print(namestr)
	foutname = today+'Cluster_DHSnormed_Module'+str(i)+'.txt'
	np.savetxt(foutname, firstcut, fmt='%i', header=theheader)

# In[ ]:



