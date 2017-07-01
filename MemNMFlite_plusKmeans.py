

# MemNMFlite_plusKmeans.py No more islands, and now we try clustering the DHSs and samples in both normalized and unnormalized iterations. 

# Since it's impossible to hold more than a few 804 x 2.5 million matrices in memory at once, i need to re-write this in a way that is perhaps computationally more efficient - done.


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

BelugaMode = True
#saveDHSmat = True


# several parameters
axis_fontsize = 30
myrandom = 0
customthresh = 0.5

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
		


AllSampleSumI_ar = []
AllDHSSumI_ar = []
bigAllDHSSum_ar = []
bigAllSampleSum_ar = []
sansvar_sums = []




errortotal = np.fabs(reconstruction - a)
error_per_sample = np.sum(errortotal, axis=1)
square_error_per_sample = np.sqrt(np.sum(np.square(errortotal), axis=1))
signal_per_sample = np.sum(a, axis=1)
recon_per_sample = np.sum(reconstruction, axis=1)

SampleStats = np.array([names,signal_per_sample,recon_per_sample,error_per_sample, square_error_per_sample]).T
foutname = today+'SampleStats_Ncomp'+str(NMFcomps)+'.txt'
np.savetxt(foutname, SampleStats, fmt="%10s %10.3f %10.3f %10.3f %10.3f")

f = open(today+'SampleCustomStats_Ncomp'+str(NMFcomps)+'.txt', 'w')
count = 0
print('starting precision recall stuff at ', mytime(), flush=True)
for sample in range(Nsamples):
	DHSar_cut = a[sample]>customthresh
	predDHSar_cut = reconstruction[sample] > customthresh
	
	
	TP = len(reconstruction[sample][DHSar_cut * predDHSar_cut]) 
	FP = len(reconstruction[sample][predDHSar_cut * np.invert( DHSar_cut)])
	TN = len(reconstruction[sample][np.invert(predDHSar_cut) * np.invert(DHSar_cut)])
	FN = len(reconstruction[sample][np.invert(predDHSar_cut) * (DHSar_cut)])

	if ((TP + FN) > 0 ):
		recall = TP / (TP + FN)
	else: 
		recall = 0
	if ((TP + FP) > 0 ):
		precision = TP / (TP + FP)
	else:
		precision=0
	accuracy = (TP + TN) /(len(reconstruction[sample]))
	print(sample, names[count], TP, FP, TN, FN, recall, precision, accuracy, file=f)	
	count +=1
f.close()


print('done with precision recall stuff ', mytime(), flush=True)
#start of outer loop
for i in range(NMFcomps):
	AllSampleSumI_ar.append([])
	AllDHSSumI_ar.append([])
	bongo = np.copy(BasisMat)
	for k in range(NMFcomps):
		if (k!=i):
			bongo[:,k]*=0
	sansvar = np.dot(bongo, MixtureMat)
	bigAllDHSSum_ar.append(np.sum(sansvar[:,0:], axis=1))
	bigAllSampleSum_ar.append(np.sum(sansvar[:,0:], axis=0))
	sansvar_sums.append(np.sum(sansvar))
	del(sansvar)
		
AllSampleSumI_ar = np.array(AllSampleSumI_ar)
AllDHSSumI_ar = np.array(AllDHSSumI_ar)
bigAllDHSSum_ar = np.array(bigAllDHSSum_ar)
bigAllSampleSum_ar = np.array(bigAllSampleSum_ar)

print('done with first loop at  ', mytime(), flush=True)

foutname = today+'samples_by_component_Ncomp'+str(NMFcomps)+'.txt'
foutnamePD = today+'samples_by_component_NcompPD'+str(NMFcomps)+'.txt'
foutnameDHSpd = today+'DHS_by_component_NcompPD'+str(NMFcomps)+'.txt'
funbigAl = pd.DataFrame(bigAllDHSSum_ar.T)
funbigSa = pd.DataFrame(bigAllSampleSum_ar.T)
funbigAl['names'] = names
np.savetxt(foutname, funbigAl.values, fmt="%10.3f "*NMFcomps+" %10s" )
funbigAl.to_csv(foutnamePD, sep='\t')
funbigSa.to_csv(foutnameDHSpd, sep='\t')


#switching to transpoe for consistency with jupyter
bigAllDHSSum_ar = bigAllDHSSum_ar.T
bigAllSampleSum_ar = bigAllSampleSum_ar.T

# time for full-scale bar graph
# unnormalized, sample-wise, adjusted via sansvar


normed_Basisums = np.sum(bigAllDHSSum_ar, axis=1)
normed_VecRats = np.copy(bigAllDHSSum_ar)
for i, rat in  enumerate(normed_VecRats):
	if (normed_Basisums[i]>0):
		rat/=normed_Basisums[i]

rawnormed_Basisums = np.sum(BasisMat, axis=1)
rawnormed_VecRats = np.copy(BasisMat)
for i, rat in  enumerate(rawnormed_VecRats):
	if (rawnormed_Basisums[i]>0):
		rat/=rawnormed_Basisums[i]



from sklearn.cluster import KMeans
# cluster unnormalized, sample-wise, adjusted via sansvar
kmeans_unnormed_sample_sansvar = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(bigAllDHSSum_ar)
normed_VecRatsS = kmeans_unnormed_sample_sansvar.cluster_centers_.T
ClusterCounts = np.bincount(kmeans_unnormed_sample_sansvar.predict(bigAllDHSSum_ar))
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_unnormed_Reweighted_Samplewise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)
print('barClusters_unnormed_Reweighted_Samplewise ',ClusterCounts )

# cluster unnormalized, sample-wise, straight from NMF
kmeans_unnormed_sample_orig = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(BasisMat)
normed_VecRatsS = kmeans_unnormed_sample_orig.cluster_centers_.T
ClusterCounts = np.bincount(kmeans_unnormed_sample_orig.predict(BasisMat)) 
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_unnormed_straight_Samplewise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)
print('barClusters_unnormed_straight_Samplewise ', ClusterCounts)




'''
#majority-based clustering ordering
WinningComponent = np.argmax(rawnormed_VecRats, axis=1)
barsortorder = np.array([])
for i in range(NMFcomps):
    barsortorder = np.append(barsortorder, np.argwhere(WinningComponent==i))
barsortorder = barsortorder.astype(int)
'''


#using k-means result to make an ordering 
WinningComponent = kmeans_unnormed_sample_sansvar.predict(bigAllDHSSum_ar)
barsortorder = np.array([])
for i in range(32):
    barsortorder = np.append(barsortorder, np.argwhere(WinningComponent==i))
barsortorder = barsortorder.astype(int)

make_bar_plot(Nrelevant=Nsamples, BarMatrix=bigAllDHSSum_ar.T, bargraph_out=today+'bars_Unnormed_Reweighted_Samplewise.png', plotClusterMode=False, myNMFcomps=NMFcomps, barsortorder=barsortorder)

make_bar_plot(Nrelevant=Nsamples, BarMatrix=normed_VecRats.T, bargraph_out=today+'bars_normed_Reweighted_Samplewise.png', plotClusterMode=False, myNMFcomps=NMFcomps, barsortorder=barsortorder)

make_bar_plot(Nrelevant=Nsamples, BarMatrix=BasisMat.T, bargraph_out=today+'bars_Unnormed_straight_Samplewise.png', plotClusterMode=False, myNMFcomps=NMFcomps, barsortorder=barsortorder)

make_bar_plot(Nrelevant=Nsamples, BarMatrix=rawnormed_VecRats.T, bargraph_out=today+'bars_normed_straight_Samplewise.png', plotClusterMode=False, myNMFcomps=NMFcomps, barsortorder=barsortorder)




# cluster normalized sample-wise straight form NMF

kmeans_normed_sample_orig = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(rawnormed_VecRats)
normed_VecRatsS = kmeans_normed_sample_orig.cluster_centers_.T
ClusterCounts = np.bincount(kmeans_normed_sample_orig.predict(rawnormed_VecRats))
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_normed_straight_Samplewise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)
print('barClusters_normed_straight_Samplewise ', ClusterCounts)


# cluster normalized sample-wise adjusted via sansvar

kmeans_normed_sample_sansvar = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(normed_VecRats)
normed_VecRatsS = kmeans_normed_sample_sansvar.cluster_centers_.T
ClusterCounts = np.bincount(kmeans_normed_sample_sansvar.predict(normed_VecRats))
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_normed_Reweighted_Samplewise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)
print('barClusters_normed_Reweighted_Samplewise ',ClusterCounts )

# use cluster ordering to do bar-plots (maybe just do that from the start)

rawnormed_Mixturesums = np.sum(MixtureMat.T, axis=1)
rawnormed_MixtureRats = np.copy(MixtureMat.T)
rawnormed_MixtureRats = rawnormed_MixtureRats[rawnormed_Mixturesums>0]
rawnormed_Mixturesums = rawnormed_Mixturesums[rawnormed_Mixturesums>0]
for i, rat in  enumerate(rawnormed_MixtureRats):
    rat/=rawnormed_Mixturesums[i]
kmeansDHS = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(rawnormed_MixtureRats)
normed_VecRatsS = kmeansDHS.cluster_centers_.T
ClusterCounts = np.bincount(kmeansDHS.predict(rawnormed_MixtureRats))
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_normed_straight_DHSwise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)
print('barClusters_normed_straight_DHSwise ', ClusterCounts)

all_kmeansDHS = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(MixtureMat.T)
normed_VecRatsS = all_kmeansDHS.cluster_centers_.T
ClusterCounts =  np.bincount(all_kmeansDHS.predict(MixtureMat.T))
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_Unnormed_straight_DHSwise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)
print('barClusters_Unnormed_straight_DHSwise ',ClusterCounts)

normed_Mixturesums = np.sum(bigAllSampleSum_ar, axis=1)
normed_MixtureRats = np.copy(bigAllSampleSum_ar)
normed_MixtureRats = normed_MixtureRats[normed_Mixturesums>0]
normed_Mixturesums = normed_Mixturesums[normed_Mixturesums>0]
for i, rat in  enumerate(normed_MixtureRats):
    rat/=normed_Mixturesums[i]

kmeansDHSw = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(bigAllSampleSum_ar)
normed_VecRatsS = kmeansDHSw.cluster_centers_.T
ClusterCounts = np.bincount(kmeansDHSw.predict(bigAllSampleSum_ar))
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_Unnormed_Reweighted_DHSwise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)

print('barClusters_Unnormed_Reweighted_DHSwise ', ClusterCounts)


kmeansDHSwnorm = KMeans(n_clusters=NMFcomps*2, random_state=0).fit(normed_MixtureRats)
normed_VecRatsS = kmeansDHSwnorm.cluster_centers_.T
ClusterCounts = np.bincount(kmeansDHSwnorm.predict(normed_MixtureRats))
make_bar_plot(Nrelevant=NMFcomps*2, BarMatrix=normed_VecRatsS, bargraph_out=today+'barClusters_normed_Reweighted_DHSwise.png', plotClusterMode=True, myNMFcomps=NMFcomps, clusterTopLabels=ClusterCounts)
print('barClusters_normed_Reweighted_DHSwise ', ClusterCounts)


