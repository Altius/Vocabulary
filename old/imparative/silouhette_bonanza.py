import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

Mixture = np.load('2017-09-25NMF_Ncomps16Mixture.npy')
cluster_assignments = np.load('2017-10-26clusterassignments36.npy')
NormedMixture =   Mixture / np.sum(Mixture, axis=0)
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from OONMFhelpers import *

today = get_today()

def define_colorsA(Ncomps, mode='newSasha'):
    if (mode=='newSasha'):
        Comp_colors = ['#A6CEE3','#1f78b4','#b2df8a','#33a02c','#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928', '#ffd700', '#AAAAAA', '#C52892', '#00bbbb']
    elif (mode=='Sasha'):
        Comp_colors = ['red', 'tan', 'lime','blue','m','k','c', 'coral', 'indigo','darkgreen','orange','grey','gold', 'lightskyblue', 'peru', 'olive']
    else:
        Comp_colors = ["#A6CEE3", "#438EC0", "#63A8A0", "#98D277", "#3BA432", "#B89B74", "#F16667", "#E62F27", "#F9A963", "#FE982C", "#ED8F47", "#C3AAD2", "#7D54A5","#B9A499", "#EAD27A" ,"#B15928"]
    if (Ncomps>16):
        np.random.seed(666)
        from matplotlib import colors as mcolors
        colornames = list(mcolors.CSS4_COLORS.keys())
        count = 16
        while (count < Ncomps):
            newcolor = colornames[np.random.randint(0,len(colornames))]
            trialcount = 0
            while ((newcolor in Comp_colors) and (trialcount < 100)):
                newcolor = colornames[np.random.randint(0,len(colornames))]
                trialcount+=1
            Comp_colors.append(newcolor)
            count+=1
    return Comp_colors
    
def make_stacked_bar_plot(Ncomps, Nrelevant, BarMatrix, bargraph_out, names = [], plotClusterMode=False, barsortorder=[], clusterTopLabels=[], colormode='newSasha', plot_title=''):
    if len(barsortorder)<1:
        barsortorder = np.arange(Nrelevant)
        #print('inventing barsortorder')
    if len(names) < 1:
        #print('inventing names')
        names = [str(i) for i in range(Nrelevant)]
        names = np.array(names)
    ttt = np.arange(Nrelevant)
    start = 0
    end = Nrelevant
    ground_pSample = ttt*0
    Comp_colors = define_colorsA(Ncomps, mode=colormode)
    plt.clf()
    plt.figure(figsize=(150,40))
    plt.bar(ttt[start:end], BarMatrix[0,start:end][barsortorder],  color=Comp_colors[0], bottom=ground_pSample[start:end], alpha=1.0)
    ground_pSample = BarMatrix[0,start:end][barsortorder]
    for i in range(1,Ncomps):
        plt.bar(ttt[start:end],BarMatrix[i,start:end][barsortorder], bottom = ground_pSample, color=Comp_colors[i], alpha=1.0)
        ground_pSample = np.sum(BarMatrix[0: i+1,start:end], axis=0)[barsortorder]
    increase_axis_fontsize()
    plt.ylabel('sum of signal in matrix',fontsize=70)
    if (len(plot_title) > 0):
        plt.title(plot_title,fontsize=70)
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
    
list_of_sils = []
for current_prototype in [j for j in range(36)]:
    prot_sils = []
    print ('doing prototype',current_prototype+1)
    cut = cluster_assignments==current_prototype
    X =  NormedMixture.T[cut]
    
    range_n_clusters = [2, 3, 4, 5, 6, 7, 8,9,10, 11, 12]
    cut = cluster_assignments==current_prototype

    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        #fig, (ax1, ax2) = plt.subplots(1, 2)
        #fig.set_size_inches(18, 7)



        # Initialize the clusterer with n_clusters value and a random generator
        # seed of 10 for reproducibility.
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        #print('done clustering', n_clusters)
        cluster_labels = (clusterer.fit(X)).labels_
        #print(cluster_labels)
        np.random.seed(666)
        desiredX = np.random.choice(X.shape[0], 5000,replace=False)
        Xcut = np.in1d(np.arange(X.shape[0]), desiredX)
        for i in range(n_clusters):
            silcut = cluster_labels[Xcut] == i
            #print('basiscut length ',len(silcut[silcut]))

            if len(silcut[silcut] ) > 40:
                chosenNum = 40
            else:
                chosenNum = len(silcut[silcut])
            possibleDHS = np.arange(len(silcut))[silcut]
            desiredDHS = np.random.choice(possibleDHS, chosenNum,replace=False)
            VisualDHScut = np.in1d(np.arange(len(silcut)), desiredDHS)
            winningcomp = np.argmax(np.sum(X[Xcut][VisualDHScut], axis=0))
            desired_order = np.argsort(-X[Xcut][VisualDHScut][:,winningcomp])    
            #desired_order = np.argsort(-NormedMixture.T[VisualDHScut][:,i])
            #UNnormedbarsortorder = get_barsortorder(Mixture.T[mixturecut])
            make_stacked_bar_plot(16, chosenNum, X[Xcut][VisualDHScut][desired_order].T, today+'p'+str(current_prototype+1)+'NormedDHSc_set_of'+str(n_clusters)+'N'+str(i+1)+'Cosine.pdf', plot_title='N'+str(i+1)+' length '+str(len(silcut[silcut])))

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        #silhouette_avg = silhouette_score(X, cluster_labels, sample_size=1000)
        silhouette_avg = silhouette_score(X[Xcut], cluster_labels[Xcut], metric='cosine')
        prot_sils.append(silhouette_avg)
        #print("For n_clusters =", n_clusters,
        #      "The average silhouette_score is :", silhouette_avg)
        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(X[Xcut], cluster_labels[Xcut], metric='cosine')

        labels_of_sils = []

        fig = plt.figure(figsize=(18,7))
        ax1=plt.gca()

        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1 but in this example all
        # lie within [-0.1, 1]
        ax1.set_xlim([-0.1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, 5000 + (n_clusters + 1) * 10])
        y_lower = 10

        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                sample_silhouette_values[cluster_labels[Xcut] == i]

            ith_cluster_silhouette_values.sort()
            labels_of_sils.append(np.median(ith_cluster_silhouette_values))

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        labels_of_sils = np.array(labels_of_sils)
        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
        plt.legend(labels_of_sils.astype(str))
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
        ax1.axvline(x=0.5, color="blue", linestyle=":")

        # 2nd Plot showing the actual clusters formed
        '''
        colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
        ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                    c=colors, edgecolor='k')

        # Labeling the clusters
        centers = clusterer.cluster_centers_
        # Draw white circles at cluster centers
        ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
                    c="white", alpha=1, s=200, edgecolor='k')

        for i, c in enumerate(centers):
            ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                        s=50, edgecolor='k')

        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("Feature space for the 1st feature")
        ax2.set_ylabel("Feature space for the 2nd feature")

        plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                      "with n_clusters = %d" % n_clusters),
                     fontsize=14, fontweight='bold')
        '''
        plt.savefig( today+'prototype'+str(current_prototype+1)+'NC'+str(n_clusters)+'silCosine.pdf')
        plt.close()
    list_of_sils.append(prot_sils)

list_of_sils = np.array(list_of_sils)
np.savetxt(today+'Cosine_list_of_sils.txt',list_of_sils)


    