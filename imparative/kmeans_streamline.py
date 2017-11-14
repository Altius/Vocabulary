import numpy as np
import pandas as pd

ClusterMode=True
if ClusterMode:
    import matplotlib
    matplotlib.use('Agg')
from OONMFhelpers import *
import matplotlib.pyplot as plt
import seaborn as sns
	
sns.set_style("whitegrid")

from sklearn.cluster import KMeans
num_clusters = 36
today = get_today()

def make_stacked_bar_plot(Ncomps, Nrelevant, BarMatrix, bargraph_out, names = [], plotClusterMode=False, barsortorder=[], clusterTopLabels=[], colormode='newSasha', plot_title=''):
    if len(barsortorder)<1:
        barsortorder = np.arange(Nrelevant)
        print('inventing barsortorder')
    if len(names) < 1:
        print('inventing names')
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
        plt.title(plot_title)
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
    

finname = '2017-09-25NMF_Ncomps16Mixture.npy'
Mixture = np.load(finname)
NormedMixture =   Mixture / np.sum(Mixture, axis=0)
kmeans_normed_DHS= KMeans(n_clusters=num_clusters, random_state=0).fit(NormedMixture.T)

clusterassignments = kmeans_normed_DHS.labels_

for i in range(num_clusters):
    print('doing ',i)
    mixturecut = clusterassignments==i
    print('basiscut length ',len(mixturecut[mixturecut]))
    
    if len(mixturecut[mixturecut] ) > 40:
        chosenNum = 40
    else:
        chosenNum = len(mixturecut[mixturecut])
    
    possibleDHS = np.arange(len(mixturecut))[mixturecut]
    desiredDHS = np.random.choice(possibleDHS, chosenNum,replace=False)
    VisualDHScut = np.in1d(np.arange(len(mixturecut)), desiredDHS)
    winningcomp = np.argmax(np.sum(NormedMixture.T[VisualDHScut], axis=0))
    desired_order = np.argsort(-NormedMixture.T[VisualDHScut][:,winningcomp])    
    #desired_order = np.argsort(-NormedMixture.T[VisualDHScut][:,i])
    #UNnormedbarsortorder = get_barsortorder(Mixture.T[mixturecut])
    make_stacked_bar_plot(16, chosenNum, NormedMixture.T[VisualDHScut][desired_order].T, today+'DHScluster'+str(i)+'outof'+str(num_clusters)+'.pdf', plot_title='N'+str(i)+' length '+str(len(mixturecut[mixturecut])))


def NMF_density_plot(mydata, foutname='', plot_title=''):
    plt.clf()
    plt.figure(figsize=(10,10))
    ax = plt.gca()
    increase_axis_fontsize()

    #ax.set_xlim(0.1, 1)
    Comp_colors = define_colorsA(16)

    currentymax = 0

    for j in range(16):
        #print('doing ',i)
        if(np.sum(NormedMixture.T[mixturecut][:,j]) > 0):
            fog = sns.distplot(NormedMixture.T[mixturecut][:,j], color=Comp_colors[j],hist=True, hist_kws={"alpha":1.0, "histtype":"step", "lw":3}, kde_kws = {"alpha":0.0})
            zag = fog.axes
            for line in zag.lines:
                yline = line.get_ydata()
                xline = line.get_xdata()
                xlinecut = xline>0.1
                if (len(yline[xlinecut]) > 0 ):
                    ymax = np.max(yline[xlinecut])
                else: 
                    ymax = 0
                if (ymax > currentymax):
                    currentymax = ymax

    currentymax*=1.1
    ax.set_xlim(0.05, 1)
    ax.set_ylim(0, currentymax)
    plt.xlabel('fractional composition')
    plt.ylabel('normalized density')
    if (len(plot_title) > 0):
        plt.title(plot_title)
    if len(foutname) > 0:
        plt.savefig(foutname)
    else:
        plt.show()
    plt.close()

np.savetxt(today+'clusterassignments'+str(num_clusters)+'.txt', clusterassignments, fmt="%i")
np.savetxt(today+'clustercenters'+str(num_clusters)+'.txt',kmeans_normed_DHS.cluster_centers_)

for i in range(num_clusters):
    mixturecut = clusterassignments == i
    NMF_density_plot(NormedMixture.T[mixturecut], foutname=today+'DensityPlot'+str(i)+'of'+str(num_clusters)+'.pdf',plot_title='N'+str(i)+' length '+str(len(mixturecut[mixturecut])))
    