ClusterMode = False

from datetime import date
import numpy as np
today = str(date.today())
if (ClusterMode):
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import pandas as pd


'''
functions:

make_significance_plot - make the signature system-wide association plots that we use to show statistical associations between components and metadata categories. 

make_significance_plot_homogeneity - make a version of these metadata associations, but test for homogeneity among NMF components 

get_rosetta - this is an important one. It figures out how to reorder the names of samples from metadata so that it matches up to the order of sample names from the masterlist matrix. Uses regular expressions. 

'''


def make_significance_plot(X, Basis, Category_str, my_rosetta, thecmap='viridis', NMFCOMPS=16, save=True, filename_addon='', PCAmode=False, write_mode=False, verbose=True, maxcats=20):
    CategoryType = X[Category_str].value_counts().keys()[0:maxcats]
    CategoryCats = X[Category_str].values[my_rosetta]

    list_of_sig = []

    for i, cat in enumerate(CategoryType[0:maxcats]):
        growthlist = [cat]
        if verbose:
            print('*****************')
            print('Category ',i,cat)
        CatCut = (CategoryCats==cat)
        growthlist+=[len(CatCut[CatCut])]
        for i in range(NMFCOMPS):
            car = mannwhitneyu(Basis[:,i][CatCut], Basis[:,i][~CatCut],alternative='greater')
            Ncats = min(maxcats, len(CategoryType[0:maxcats]))
            adjustedp = (car[1]+1e-30)*NMFCOMPS*Ncats
            growthlist += [-1*np.log10(adjustedp)]
        list_of_sig.append(growthlist)
        
    colnames = [Category_str.replace(" ", ""), 'Count'] + ['Comp'+str(i+1) for i in range(NMFCOMPS)]
    
    CategoryChart = pd.DataFrame(list_of_sig, columns=colnames)
    if write_mode:
        CategoryChart.to_csv(filename_addon +  Category_str+'MWmatrix.csv', sep='\t') 
    CategoryChartMatrix = CategoryChart.values[:,2:].astype(float)

    plt.clf()
    
    myfs = 45
    plt.figure(figsize=(35,len(CategoryType[0:maxcats])*2))
    plt.imshow(CategoryChartMatrix, cmap=thecmap, vmin=-3, vmax=30)
    if (PCAmode):
        plt.xlabel('Principal Component',fontsize=myfs)
    else:
        plt.xlabel('NMF component',fontsize=myfs)
    plt.ylabel(Category_str,fontsize=myfs)
    plt.yticks(np.arange(Ncats), CategoryType[0:maxcats], rotation='horizontal',fontsize=myfs)
    plt.xticks(np.arange(NMFCOMPS), (np.arange(NMFCOMPS)+1).astype(str), rotation='vertical',fontsize=myfs)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label(r'- $\log_{10} (p*$'+str(NMFCOMPS)+r'$*$'+str(Ncats)+r'$)$',fontsize=myfs)
    cbar_ax = cbar.ax
    cbar_ax.tick_params(labelsize=myfs)
    for i in cbar_ax.get_yticklabels():
        i.set_fontsize(myfs)
    if(save):
        plt.savefig(filename_addon +  Category_str+'MWplot.pdf', bbox_inches='tight')
    plt.show()

    return (CategoryChartMatrix, CategoryType[0:maxcats])
    

def make_significance_plot_WSO(X, Basis, Category_str, my_rosetta, thecmap='binary', NMFCOMPS=16, save=True, filename_addon='', PCAmode=False, write_mode=False, verbose=True, maxcats=20):
    CategoryType = X[Category_str].value_counts().keys()[0:maxcats]
    CategoryCats = X[Category_str].values[my_rosetta]
    
    CategoryType = np.sort(CategoryType)

    list_of_sig = []

    for i, cat in enumerate(CategoryType[0:maxcats]):
        growthlist = [cat]
        if verbose:
            print('*****************')
            print('Category ',i,cat)
        CatCut = (CategoryCats==cat)
        growthlist+=[len(CatCut[CatCut])]
        for i in range(NMFCOMPS):
            car = mannwhitneyu(Basis[:,i][CatCut], Basis[:,i][~CatCut],alternative='greater')
            Ncats = min(maxcats, len(CategoryType[0:maxcats]))
            adjustedp = (car[1]+1e-30)*NMFCOMPS*Ncats
            growthlist += [-1*np.log10(adjustedp)]
        list_of_sig.append(growthlist)
        
    colnames = [Category_str.replace(" ", ""), 'Count'] + ['Comp'+str(i+1) for i in range(NMFCOMPS)]
    
    CategoryChart = pd.DataFrame(list_of_sig, columns=colnames)
    if write_mode:
        CategoryChart.to_csv(filename_addon +  Category_str+'MWmatrix.csv', sep='\t') 
    CategoryChartMatrix = CategoryChart.values[:,2:].astype(float)

    plt.clf()
    
    myfs = 45
    plt.figure(figsize=(35,len(CategoryType[0:maxcats])*2))
    plt.imshow(CategoryChartMatrix, cmap=thecmap, vmin=-3, vmax=30)
    if (PCAmode):
        plt.xlabel('Principal Component',fontsize=myfs)
    else:
        plt.xlabel('NMF component',fontsize=myfs)
    plt.ylabel(Category_str,fontsize=myfs)
    plt.yticks(np.arange(Ncats), CategoryType[0:maxcats], rotation='horizontal',fontsize=myfs)
    plt.xticks(np.arange(NMFCOMPS), (np.arange(NMFCOMPS)+1).astype(str), rotation='vertical',fontsize=myfs)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label(r'- $\log_{10} (p*$'+str(NMFCOMPS)+r'$*$'+str(Ncats)+r'$)$',fontsize=myfs)
    cbar_ax = cbar.ax
    cbar_ax.tick_params(labelsize=myfs)
    for i in cbar_ax.get_yticklabels():
        i.set_fontsize(myfs)
    if(save):
        plt.savefig(filename_addon +  Category_str+'MWplot.pdf', bbox_inches='tight')
    plt.show()

    return (CategoryChartMatrix, CategoryType[0:maxcats])
        


def make_significance_plot_homogeneity(X, homogeneity, Category_str, my_rosetta, thecmap='viridis', NMFCOMPS=16, save=True,filename_addon='', verbose=True, maxcats=20):
    CategoryType = X[Category_str].value_counts().keys()[0:maxcats]
    CategoryCats = X[Category_str].values[my_rosetta]

    list_of_sig = []

    for i, cat in enumerate(CategoryType[0:maxcats]):
        growthlist = [cat]
        if verbose:
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
        if verbose:
            print ('mean homogeneity of ',cat,np.mean(homogeneity[CatCut]))
            print ('mean homogeneity of anti-',cat,np.mean(homogeneity[~CatCut]))
        Ncats = min(maxcats, len(CategoryType[0:maxcats]))
        adjustedp = (car[1]+1e-30)*1*Ncats
        growthlist += [-1*np.log10(adjustedp)]




        list_of_sig.append(growthlist)
        
    colnames = [Category_str.replace(" ", ""), 'Count'] +['A']
    
    CategoryChart = pd.DataFrame(list_of_sig, columns=colnames)
    CategoryChartMatrix = CategoryChart.values[:,2:].astype(float)

    plt.clf()
    plt.figure(figsize=(len(CategoryType[0:maxcats])*2,6))
    plt.imshow(CategoryChartMatrix.T, cmap=thecmap, vmin=-3, vmax=15)

    plt.xlabel(Category_str,fontsize=25)
    plt.xticks(np.arange(Ncats), CategoryType[0:maxcats], rotation='vertical',fontsize=25)

    cbar = plt.colorbar(fraction=0.046, pad=0.04, ticklabel_size=24)
    cbar.set_label(r'- $\log_{10} (p*$'+str(Ncats)+r'$)$',fontsize=25)
    cbar_ax = cbar.ax
    cbar_ax.tick_params(labelsize=35)
    for i in cbar_ax.get_yticklabels():
        i.set_fontsize(35)

    if(save):
        plt.savefig(filename_addon +  Category_str+'MWhom_plot.pdf')
    plt.show()

    return (CategoryChartMatrix, CategoryType[0:maxcats])


def get_rosetta(MetaDataMat, names):
    import re
    DSnos = []
    DSnos_naked = []
    for name in names:
        p = re.search('(.)(DS\w+)', name)
        DSnos.append(p.group(2))
        DSnos_naked.append(p.group(2)[2:])
    rosetta = []
    for i in DSnos:
        rosetta.append(np.argwhere(i == MetaDataMat['DS_plus'].values)[0][0])
    return rosetta


