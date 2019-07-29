ClusterMode = False

import numpy as np
import pandas as pd
if (ClusterMode):
	import matplotlib
	matplotlib.use('Agg')
	
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec




def return_expected_matrix(Gcut, maxcolor, secondmaxcolor, secondmaxcolor_vals, strings_of_labels,  rownormed=False):

    bincount_matrix_cutoff = []
    for comp in range(16):
        bincut_row = []
        maincut = maxcolor[Gcut] == comp
        for j in range(16):
            seccut = (secondmaxcolor[Gcut] == j) *  (secondmaxcolor_vals[Gcut]>0)
            if comp == j:
                seccut = secondmaxcolor_vals[Gcut]==0
            N_in_square_cut = maincut * seccut
            N_in_square = np.sum(N_in_square_cut)
            bincut_row.append(N_in_square)
        bincount_matrix_cutoff.append(bincut_row)
    bincount_matrix_cutoff = np.array(bincount_matrix_cutoff)
    if rownormed:
        bincount_matrix_cutoff_rownormed = (bincount_matrix_cutoff.T / np.sum(bincount_matrix_cutoff, axis=1)).T
        bincount_matrix_cutoff_rownormedPD = pd.DataFrame(bincount_matrix_cutoff_rownormed, columns=strings_of_labels)
        bincount_matrix_cutoff_rownormedPD['domcomp'] = strings_of_labels
        bincount_matrix_cutoff_rownormedPD.set_index('domcomp', inplace=True)
    
    bincount_matrix_cutoffPD = pd.DataFrame(bincount_matrix_cutoff, columns=strings_of_labels)
    bincount_matrix_cutoffPD['domcomp'] = strings_of_labels
    bincount_matrix_cutoffPD.set_index('domcomp', inplace=True)

        #return both  counts and row-normed version
    if rownormed:
        return [bincount_matrix_cutoffPD, bincount_matrix_cutoff_rownormedPD]

    return bincount_matrix_cutoffPD

# not the prettiest function, but 
def big_grid_plot(Gcut, maxcolor, secondmaxcolor, maxcolor_vals, secondmaxcolor_vals, WSO, color_scheme, strings_of_labels, plt_title = '', do_row_normalized_shading=True, do_extra_bargraphs=True, logOOEmode=False, adjust_significance=False, ExpectedMatrix = [], pmatrix = [] ):
    maxcolor_bincount = np.bincount(maxcolor[Gcut], minlength=16)
    seccolor_bincount = np.bincount(secondmaxcolor[Gcut], minlength=16)

    seccut = secondmaxcolor_vals>0
    # possibly need to revise this - by default, the expectation value is just all sites
    maxcolor_bincount_OOE = np.bincount(maxcolor[Gcut], minlength=16) /  np.bincount(maxcolor, minlength=16) * (len(maxcolor) / len(maxcolor[Gcut]))

    seccolor_bincount_OOE =  np.bincount(secondmaxcolor[Gcut*seccut], minlength=16) /  np.bincount(secondmaxcolor[seccut], minlength=16) * (len(secondmaxcolor[seccut]) / len(secondmaxcolor[Gcut* seccut]))
    
    
    if do_row_normalized_shading:
        [bincount_matrix_cutoff_G_PD, bincount_matrix_cutoff_G_rownormedPD] = return_expected_matrix(Gcut,maxcolor, secondmaxcolor, secondmaxcolor_vals, strings_of_labels, rownormed=True)
    else:
        bincount_matrix_cutoff_G_PD = return_expected_matrix(Gcut,maxcolor, secondmaxcolor,secondmaxcolor_vals, strings_of_labels)
    

    
    gs1 = gridspec.GridSpec(5,5)
    fig = plt.figure(figsize = (24,22))
    plt.clf()
    f = plt.subplot2grid((5,5), (0,1), colspan=3, rowspan=3)
    
    
    if (do_row_normalized_shading):
        kox= sns.heatmap(bincount_matrix_cutoff_G_rownormedPD[strings_of_labels[WSO]].reindex(strings_of_labels[WSO]), annot=bincount_matrix_cutoff_G_PD[strings_of_labels[WSO]].reindex(strings_of_labels[WSO]).values.astype(int), cmap='binary', vmax=0.5, fmt='d', cbar_kws={'fraction':0.046, 'pad':0.04} )
        kox.collections[0].colorbar.set_label('secondary fraction', fontsize=25)
    elif logOOEmode:
        if len(ExpectedMatrix) < 1: 
            print('error, must pass return matrix')
            plt.close()
            return
        group_rat = bincount_matrix_cutoff_G_PD / ExpectedMatrix
        group_fix = ExpectedMatrix.sum().sum() / bincount_matrix_cutoff_G_PD.sum().sum()

        thing_to_heatmap = np.log10(group_fix*group_rat[strings_of_labels[WSO]].reindex(strings_of_labels[WSO])+1e-5) / np.log10(2) #using log2 
        
        if (adjust_significance):            
            temp_annotations = thing_to_heatmap.values
            annotations = []
            
            for i in range(len(temp_annotations)):
                irow = []
                for j in range(len(temp_annotations)):
                    cellvalue = '{:6.2f}'.format(temp_annotations[i][j]).strip(' ')
                    if pmatrix[i][j]:
                        irow.append(cellvalue+'*')
                    else:
                        irow.append(cellvalue)
                annotations.append(irow)
            annotations = np.array(annotations)
            
            #annotations = temp_annotations
            annotationsPD = pd.DataFrame(annotations, index =strings_of_labels[WSO], columns = strings_of_labels[WSO] )
            print(annotations.shape)
            kox= sns.heatmap(thing_to_heatmap, annot=annotationsPD, cmap='RdYlGn', vmax=8, vmin=-8, fmt='s', cbar_kws={'fraction':0.046, 'pad':0.04}, annot_kws={'fontsize':12}, mask =  thing_to_heatmap<-10 )
            
        else:
            kox= sns.heatmap(thing_to_heatmap, annot=True, cmap='RdYlGn', vmax=8, vmin=-8, fmt='6.2f', cbar_kws={'fraction':0.046, 'pad':0.04}, mask =  thing_to_heatmap<-10 )


        kox.collections[0].colorbar.set_label('log2 O/E', fontsize=25)

    
    
    ax =plt.gca()
    kox.collections[0].colorbar.ax.tick_params(labelsize=25)
    kox.set_aspect("equal")
    the_fontsize=24
    ax.set_ylabel('primary component',fontsize=20)    
    ax.set_xlabel('secondary component', fontsize=20)
    ax.get_xticklabels()
    for n, i in enumerate(ax.get_xticklabels()):
        i.set_color(color_scheme[n])
        i.set_fontsize(the_fontsize)
    for n, i in enumerate(ax.get_yticklabels()):
        i.set_fontsize(the_fontsize)
        i.set_color(color_scheme[n])
    plt.title(plt_title, fontsize=40)
    

    if (do_extra_bargraphs):
        
        #left plot for O/E for main comp 
        
        leftplot = plt.subplot2grid((5,5), (0,0), rowspan=3, position=[-0.05, 0.395, 0.1, 0.525])

        leftplot.barh(['C'+str(i+1) for i in range(16)], maxcolor_bincount_OOE[WSO][::-1], color=color_scheme[::-1], tick_label=strings_of_labels[WSO][::-1])
        leftplot.plot([1, 1], [-1, 16], '--k')
        leftplot.plot([2, 2], [-1, 16], '--k')

        leftplot_ax = plt.gca()
        leftplot_ax.set_xticks([1 , 2])
        leftplot_ax.set_yticklabels([])
        leftplot_ax.set_xticklabels(['1', '2'])
        leftplot_ax.get_yaxis().set_visible(False)
        for n, i in enumerate(leftplot_ax.get_xticklabels()):
            i.set_fontsize(20)
        leftplot_ax.invert_xaxis()
        plt.xlabel('O/E', fontsize=24)
        plt.box(False)

        #bottom plot for secondary O/E
        botplot = plt.subplot2grid((5,5), (4,0), colspan=3, position=[0.255, 0.1, 0.485, 0.1])
        botplot.bar(['C'+str(i+1) for i in range(16)], seccolor_bincount_OOE[WSO], color=color_scheme, tick_label=strings_of_labels[WSO])
        botplot.plot( [-1, 16], [1, 1], '--k')
        botplot.plot( [-1, 16], [2, 2], '--k')
        botplot_ax = plt.gca()
        botplot_ax.set_yticks([1 , 2])
        botplot_ax.set_xticklabels([])
        botplot_ax.set_yticklabels(['1', '2'])
        botplot_ax.get_xaxis().set_visible(False)
        for n, i in enumerate(botplot_ax.get_yticklabels()):
            i.set_fontsize(20)
        botplot_ax.invert_yaxis()
        plt.ylabel('O/E', fontsize=24)
        plt.box(False)
        
    plt.savefig('test.pdf', bbox_inches='tight', transparent=True)

    plt.show()

    plt.close()