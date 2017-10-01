ClusterMode = True

from datetime import date
import numpy as np
today = str(date.today())
if (ClusterMode):
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

import pandas as pd

def make_significance_plot(X, Basis, Category_str, my_rosetta, thecmap='viridis', NMFCOMPS=16, save=True, filename_addon=''):
    CategoryType = X[Category_str].value_counts().keys()[0:20]
    CategoryCats = X[Category_str].values[my_rosetta]

    list_of_sig = []

    for i, cat in enumerate(CategoryType[0:20]):
        growthlist = [cat]
        print('*****************')
        print('Category ',i,cat)
        CatCut = (CategoryCats==cat)
        growthlist+=[len(CatCut[CatCut])]
        for i in range(NMFCOMPS):
            car = mannwhitneyu(Basis[:,i][CatCut], Basis[:,i][~CatCut],alternative='greater')
            #print(car)
            Ncats = min(20, len(CategoryType[0:20]))
            adjustedp = (car[1]+1e-30)*NMFCOMPS*Ncats
            growthlist += [-1*np.log10(adjustedp)]
        list_of_sig.append(growthlist)
        
    colnames = [Category_str.replace(" ", ""), 'Count'] + ['Comp'+str(i) for i in range(NMFCOMPS)]
    
    CategoryChart = pd.DataFrame(list_of_sig, columns=colnames)
    CategoryChartMatrix = CategoryChart.values[:,2:].astype(float)

    plt.clf()
    plt.figure(figsize=(30,len(CategoryType[0:20])*2))
    plt.imshow(CategoryChartMatrix, cmap=thecmap, vmin=-3, vmax=30)
    plt.xlabel('NMF component',fontsize=25)
    plt.ylabel(Category_str,fontsize=25)
    plt.yticks(np.arange(Ncats), CategoryType[0:20], rotation='horizontal',fontsize=25)
    plt.xticks(np.arange(NMFCOMPS), np.arange(NMFCOMPS).astype(str), rotation='vertical',fontsize=25)
    cbar = plt.colorbar()
    cbar.set_label(r'- $\log_{10} (p*$'+str(NMFCOMPS)+r'$*$'+str(Ncats)+r'$)$',fontsize=25)
    if(save):
        plt.savefig(today+Category_str+'MWplotNC'+str(NMFCOMPS)+filename_addon+'.pdf')
    else:
        plt.show()

    return (CategoryChartMatrix, CategoryType[0:20])

def make_significance_plot_homogeneity(X, homogeneity, Category_str, my_rosetta, thecmap='viridis', NMFCOMPS=16, save=True,filename_addon=''):
    CategoryType = X[Category_str].value_counts().keys()[0:20]
    CategoryCats = X[Category_str].values[my_rosetta]

    list_of_sig = []

    for i, cat in enumerate(CategoryType[0:20]):
        growthlist = [cat]
        print('*****************')
        print('Category ',i,cat)
        CatCut = (CategoryCats==cat)
        growthlist+=[len(CatCut[CatCut])]
        majorcount = 0
        semimajorcount = 0
        bigcount=0
        minorcount = 0
        noncount = 0
        car = mannwhitneyu(homogeneity[CatCut], homogeneity[~CatCut],alternative='greater')
        print ('mean homogeneity of ',cat,np.mean(homogeneity[CatCut]))
        print ('mean homogeneity of anti-',cat,np.mean(homogeneity[~CatCut]))
        #print(car)
        Ncats = min(20, len(CategoryType[0:20]))
        adjustedp = (car[1]+1e-30)*1*Ncats
        growthlist += [-1*np.log10(adjustedp)]




        list_of_sig.append(growthlist)
        
    colnames = [Category_str.replace(" ", ""), 'Count'] +['A']
    
    CategoryChart = pd.DataFrame(list_of_sig, columns=colnames)
    CategoryChartMatrix = CategoryChart.values[:,2:].astype(float)

    plt.clf()
    plt.figure(figsize=(len(CategoryType[0:20])*2,6))
    plt.imshow(CategoryChartMatrix.T, cmap=thecmap, vmin=-3, vmax=15)
    #plt.xlabel('NMF component',fontsize=25)
    plt.xlabel(Category_str,fontsize=25)
    plt.xticks(np.arange(Ncats), CategoryType[0:20], rotation='vertical',fontsize=25)
    #plt.xticks(np.arange(NMFCOMPS), np.arange(NMFCOMPS).astype(str), rotation='vertical',fontsize=25)
    cbar = plt.colorbar()
    cbar.set_label(r'- $\log_{10} (p*$'+str(Ncats)+r'$)$',fontsize=25)
    if(save):
        plt.savefig(today+Category_str+'homogMWplotNC'+str(NMFCOMPS)+filename_addon+'.pdf')
    else:
        plt.show()

    return (CategoryChartMatrix, CategoryType[0:20])
#come back here

def get_rosetta(MetaDataMat, names):
    import re
    DSnos = []
    DSnos_naked = []
    for name in names:
        p = re.search('(_)(DS\w+)', name)
        DSnos.append(p.group(2))
        DSnos_naked.append(p.group(2)[2:])
    rosetta = []
    for i in DSnos:
        rosetta.append(np.argwhere(i == MetaDataMat['DS_plus'].values)[0][0])
    return rosetta


