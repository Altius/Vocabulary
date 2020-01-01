# Alexander Murartov
# Altius Insstitute 
# August 2019

# this code assigns regulatory domains to genes according to a set of parameters and evaluates their "uniqueness" and "completeness" ratios. 
# The values of the parameters are read from a CSV file, with the row number determined by a command line argument

import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import sys


testmode=False
write_output=True
today = '81519_coarse_screen'

if len(sys.argv) < 2:
    print('syntax blah.py N')
    sys.exit()

line_no = int(sys.argv[1])
#read paramter files 
table = pd.read_csv('80519_parameter_list_linear.csv', sep='\t')
params = table.iloc[line_no]
thef1 = params.f1
thef2 = params.f2
thef3 = params.f3
thef4 = params.f4
thef5 = params.f5
thecosthresh = params.costhresh
print('parameters ',thef1,thef2,thef3,thef4,thef5,thecosthresh)

#parameter set 1

#thef1=1; thef2=1; thef3=1; thef4=1; thef5=2;

#parameter set 2
#thef1 = 2*np.sqrt(4/7); thef2 = np.sqrt(4/7); thef3=np.sqrt(4/7); thef4=np.sqrt(4/7); thef5=2;

#parameter set 3
#thef1 = np.sqrt(4/7); thef2 = 2*np.sqrt(4/7); thef3=np.sqrt(4/7); thef4=np.sqrt(4/7); thef5=2;

chr_list = [str(i) for i in range(1,23)] +['X', 'Y']
if testmode:
    chr_list=['16']

novanilla = True


# read pre-computed gene csv file
geneDF = pd.read_csv('80419_ensemble93_genes.csv', sep='\t', index_col=0)
TSS = np.copy(geneDF.end.values)
strand_plus_cut = geneDF.strand.values=='+'
TSS[strand_plus_cut] = geneDF.start.values[strand_plus_cut]
geneDF['TSS'] = TSS

ML = pd.read_csv('masterlist_DHSs_733samples_WM20180608_all_indexIDs.txt', header=None, names=['chr', 'start', 'end', 'id', 'ML4', 'ML5', 'ML6', 'ML7', 'ML8', 'ML9', 'ML10'], dtype = str, sep='\t')

DHSweights = ML.ML4.values.astype(float) /ML.ML5.values.astype(float)
Mixture = np.load('2018-06-08NC16_NNDSVD_Mixture.npy').T
NormedMixture = (Mixture.T / np.sum(Mixture, axis=1)).T
ML['middle'] = (ML.start.values.astype(int) + ML.end.values.astype(int))/2

# my weighted KDE function.
def wegihted_kde(start, end, peak_centers, peak_bandwidths, weights, comp_fractions, extra_distance=0, nbins=1000):
    X_plot = np.linspace(start-extra_distance, end+extra_distance, nbins)
    my_dens = np.zeros(len(X_plot))
    for j, DHScord in enumerate(peak_centers):    
        microdens = np.exp(-1* np.square(DHScord - X_plot) / (2*peak_bandwidths[j]**2)  )
        my_dens += microdens* weights[j]* comp_fractions[j]
    return ([X_plot, my_dens])
    
#     
def get_chrom_features(the_chr='chrX'):
    short_chr = the_chr.lstrip('chr')
    gene_min_x = 0
    gene_max_x = ML[ML.chr==the_chr]['end'].values.astype(int).max() //1e6 + 2
    ML_gene_cut = (ML.chr.values==the_chr)* (ML.end.values.astype(float) > gene_min_x*1e6) * (ML.start.values.astype(float) < gene_max_x*1e6)
    gene_Mixture = NormedMixture[ML_gene_cut]
    DHS_X = ML['middle'].values.astype(float)[ML_gene_cut]/1e6
    return [DHS_X, gene_Mixture, gene_max_x]
# compute likelihood     
def compute_likelihood_bodyDHS(comb_geneDF, the_gene ='GATA1', dist_left = 1, dist_right=1, costhresh=0.5, f1=1, f2=1, f3=1, f4=1, f5=2, noplot=False, DHS_X=[], gene_Mixture=[], gene_max_x=[], compute_on_fly=True, verbose=False, addon='', full_return=False):
    the_start = comb_geneDF.loc[the_gene].start/1e6 - dist_left/1000
    the_end = comb_geneDF.loc[the_gene].end/1e6 + dist_right/1000
    the_chr = 'chr'+comb_geneDF.loc[the_gene]['chr']
    if compute_on_fly:
        DHS_X, gene_Mixture, gene_max_x = get_chrom_features(the_chr)  
    
    cut = (DHS_X > the_start) * (DHS_X < the_end) 
    if verbose:
        print(len(cut[cut]))
    
    temp_Mixture = gene_Mixture[cut]
    
    if novanilla:
        L=0
    else:    
        the_size = the_end - the_start  
        gene_size = comb_geneDF.loc[the_gene].end/1e6 - comb_geneDF.loc[the_gene].start/1e6
        the_incoherence = np.sum(np.std(temp_Mixture, axis=0))/ np.sqrt(15)
        the_coherence = 1-the_incoherence
        foreground_cut = DHS_X[cut] > 0
        N_foreground = len(foreground_cut[foreground_cut])
        if (N_foreground<1):
            f_foreground=0
            true_start = the_start
            true_end = the_end
        else:
            f_foreground = len(foreground_cut[foreground_cut])/len(foreground_cut)
            true_start = np.min(DHS_X[cut][foreground_cut])
            true_end = np.max(DHS_X[cut][foreground_cut])

        center_of_mass = np.mean(DHS_X[cut][foreground_cut])
        gene_body_cut = (DHS_X[cut] > comb_geneDF.loc[the_gene].start/1e6) * (DHS_X[cut] < comb_geneDF.loc[the_gene].end/1e6) 
        body_DHS = max([len(gene_body_cut[foreground_cut*gene_body_cut]), 1])

    

        if (center_of_mass <true_start):
            d_com = (comb_geneDF.loc[the_gene].start/1e6 - center_of_mass) / ( comb_geneDF.loc[the_gene].start/1e6 - true_start)
        elif (center_of_mass > comb_geneDF.loc[the_gene].end/1e6):
            d_com = (center_of_mass - comb_geneDF.loc[the_gene].end/1e6)/(true_end - comb_geneDF.loc[the_gene].end/1e6)
        else:
            d_com = 0
        if d_com < 0:
            d_com = 0
        
        W1 =  min([1, (N_foreground/(f5*body_DHS))])
        W2 = the_coherence
        W3 = f_foreground
        W4 = (1-d_com)
        if verbose:
            print('weight of f1 (size)',W1, 'f1=',f1)
            print('weight of f2 (coherence)',W2, 'f2=',f2)
            print('weight of f3 (f_foreground)',W3, 'f3=',f3)
            print('weight of f4 (centeredness)',W4, 'f4=',f4)

        L = f1*W1 + f2*W2 + f3*W3 + f4*W4
    my_extra = 10 + 0.05*(the_end-the_start)*1000
    if not noplot:        
        plot_gene_miniregion_with_extra(DHS_X,  gene_max_x, gene_Mixture, the_start, the_end, the_chr=the_chr, foutname=today+'_closeup'+the_gene+'_'+addon+'.pdf', showplot=True, region_gene = the_gene, extra_d=my_extra )
    

    
    #perturbed L 
    foreground_cut = cdist(temp_Mixture, [np.mean (temp_Mixture, axis=0)], metric='cosine') < costhresh
    foreground_cut = foreground_cut.flatten()
    #recalculate size? maybe but i won't do it for now
    N_foreground = len(foreground_cut[foreground_cut])
    if (N_foreground<1):
        f_foreground=0
    else:
        f_foreground = len(foreground_cut[foreground_cut])/len(foreground_cut)
    gene_body_cut = (DHS_X[cut] > comb_geneDF.loc[the_gene].start/1e6) * (DHS_X[cut] < comb_geneDF.loc[the_gene].end/1e6) 
    
    body_DHS = max([len(gene_body_cut[foreground_cut*gene_body_cut]), 1])
    if (N_foreground<1):
        the_coherence = 0
        true_start =the_start
        true_end = the_end

    else:
        the_incoherence = np.sum(np.std(temp_Mixture[foreground_cut], axis=0))/ np.sqrt(15)
        the_coherence = 1-the_incoherence
        true_start = np.min(DHS_X[cut][foreground_cut])
        true_end = np.max(DHS_X[cut][foreground_cut])

    center_of_mass = np.mean(DHS_X[cut][foreground_cut])


    if (center_of_mass <true_start):
        d_com = (comb_geneDF.loc[the_gene].start/1e6 - center_of_mass) / ( comb_geneDF.loc[the_gene].start/1e6 - true_start)
    elif (center_of_mass > comb_geneDF.loc[the_gene].end/1e6):
        d_com = (center_of_mass - comb_geneDF.loc[the_gene].end/1e6)/(true_end - comb_geneDF.loc[the_gene].end/1e6)
    else:
        d_com = 0
    if d_com < 0:
        d_com = 0
        
    W1 =  min([1, (N_foreground/(f5*body_DHS))])
    W2 = the_coherence
    W3 = f_foreground
    W4 = (1-d_com)
    if verbose:
        print('weight of f1 (size)',W1, 'f1=',f1)
        print('weight of f2 (coherence)',W2, 'f2=',f2)
        print('weight of f3 (f_foreground)',W3, 'f3=',f3)
        print('weight of f4 (centeredness)',W4, 'f4=',f4)

    perturbedL =  f1*W1 + f2*W2 + f3*W3 + f4*W4
    if not noplot:        
        plot_gene_miniregion_with_extra(DHS_X,  gene_max_x, gene_Mixture, the_start, the_end, the_chr=the_chr, foutname=today+'_closeup_perturbed'+the_gene+'_'+addon+'.pdf', showplot=True, region_gene = the_gene, the_foreground=foreground_cut, extra_d=my_extra )


    if full_return:
        return_dict = {}
        return_dict['gene'] = the_gene
        return_dict['chrom'] = the_chr
        return_dict['vanilla_likelihood'] = L
        return_dict['final_likelihood'] = perturbedL
        return_dict['final_mixture'] = np.mean(gene_Mixture[cut][foreground_cut], axis=0)
        return_dict['final_DHS_indices_owned'] = np.arange(len(DHS_X))[cut][foreground_cut]
        return_dict['final_DHS_indices_rejected'] = np.arange(len(DHS_X))[cut][~foreground_cut]
        return_dict['final_probs'] = np.array([W1, W2, W3, W4])
        return return_dict
    else:
        return [L, perturbedL, N_foreground]

# search for best likelihood based on a pre-set list of possible extensions in either direction        
def search_best_likelihood(comb_geneDF, gene ='MYH7', addon='', noplot=False,DHS_X=[], gene_Mixture=[], gene_max_x=[], compute_on_fly=True, f1=1, f2=1, f3=1, f4=1, f5=2,costhresh=0.5, search_space= [0,1,2,5,10,20,30,40,50,60,70,80,90,100]):
    the_chr = 'chr'+comb_geneDF.loc[gene]['chr']
    #print('searching ',gene)
    if compute_on_fly:
        DHS_X, gene_Mixture, gene_max_x = get_chrom_features(the_chr)  
    list_of_iterations = []
    for i in search_space:
        for j in search_space:
            results = compute_likelihood_bodyDHS(comb_geneDF, noplot=True, dist_left=i, dist_right=j, the_gene=gene, f1=f1, f2=f2, f3=f3,  f4=f4, f5=f5, costhresh=costhresh, DHS_X=DHS_X, gene_Mixture=gene_Mixture, gene_max_x=gene_max_x, compute_on_fly=False)
            list_of_iterations.append([i,j]+results)
    list_of_iterationsPD = pd.DataFrame(list_of_iterations, columns=['left_kb', 'right_kb', 'init_likelihood', 'adjusted_likelihood', 'N_foreground'])
    
    best_iteration = list_of_iterationsPD.loc[list_of_iterationsPD.adjusted_likelihood.idxmax()].copy()
    best_likelihood =  best_iteration.adjusted_likelihood
    bestcut = list_of_iterationsPD.adjusted_likelihood.values == best_likelihood
    if np.sum(bestcut) > 1:
        list_of_iterationsPD['total_extension'] = list_of_iterationsPD.left_kb + list_of_iterationsPD.right_kb
        best_foreground = list_of_iterationsPD[bestcut].N_foreground.max()
        mostcut = list_of_iterationsPD[bestcut].N_foreground == best_foreground
        #print(list_of_iterationsPD[bestcut].values)
        if np.sum(mostcut) > 1:       
            best_iteration = list_of_iterationsPD[bestcut][mostcut].loc[list_of_iterationsPD[bestcut][mostcut].total_extension.idxmin()].copy()
        else:
            best_iteration = list_of_iterationsPD[bestcut].loc[list_of_iterationsPD[bestcut].total_extension.idxmax()].copy()
        
    return list_of_iterationsPD, best_iteration


total_assigned_DHSs = []
unique_assigned_DHSs = []
total_DHS_per_chrom = []

# perform procedure for each chromosome
for chrom in chr_list:
    geneDFcut = ((geneDF['biotype']=='protein_coding').values * (geneDF['chr']== chrom).values)
    prev_chr = ''
    best_domain_properties_table_SET1 = []
    best_domain_properties_table_SET1_table = []
    
    gene_DHS_ownership_table = []
    gene_DHS_rejection_table = []
    gene_mixture_table = []
    gene_performance_table = []
    
    
    for i in range(geneDF[geneDFcut].shape[0]):
        gene = geneDF[geneDFcut].iloc[i].name
        print('doing ',gene)
        #the_topcomp = geneDF[geneDFcut].iloc[i].top_comp.astype(int).astype(str)
        the_chr = geneDF[geneDFcut].iloc[i]['chr']
        if the_chr!=prev_chr:
            print('switching chr', the_chr)
            DHS_X, gene_Mixture, gene_max_x = get_chrom_features('chr'+the_chr)
            prev_chr = the_chr
        likelihoods, best_iteration  = search_best_likelihood(geneDF, gene, noplot=True, DHS_X=DHS_X, gene_Mixture=gene_Mixture, gene_max_x=gene_max_x, compute_on_fly=False, f1=thef1, f2=thef2, f3=thef3,f4=thef4,f5=thef5,costhresh=thecosthresh)
        
        if testmode:
            print(likelihoods.sort_values(by='adjusted_likelihood', ascending=False).head())
            
        
        best_left_kb = best_iteration.left_kb
        best_right_kb = best_iteration.right_kb
        
        
        
        best_domain_properties = compute_likelihood_bodyDHS(geneDF, noplot=True, dist_left=best_left_kb, dist_right=best_right_kb, the_gene=gene, f1=thef1, f2=thef2, f3=thef3,f4=thef4,f5=thef5,costhresh=thecosthresh, compute_on_fly=False,  full_return=True, DHS_X=DHS_X, gene_Mixture=gene_Mixture, gene_max_x=gene_max_x)
        best_domain_properties['left_kb'] = best_left_kb
        best_domain_properties['right_kb'] = best_right_kb
        
        gene_DHS_indices = best_domain_properties['final_DHS_indices_owned']
        for DHSnum in gene_DHS_indices:
            gene_DHS_ownership_table.append([gene, DHSnum])
            
        gene_DHS_indices_rejected = best_domain_properties['final_DHS_indices_rejected']
        for DHSnum in gene_DHS_indices_rejected:
            gene_DHS_rejection_table.append([gene, DHSnum])
            
        gene_mixture_table.append([gene] + best_domain_properties['final_mixture'].tolist())
        gene_performance_table.append([gene] + best_domain_properties['final_probs'].tolist())
        
        gene_size_kbp = geneDF.loc[gene].end/1e3 - geneDF.loc[gene].start/1e3
        best_domain_properties_table_SET1_table.append([gene,gene_size_kbp, best_left_kb,best_right_kb,likelihoods.loc[likelihoods.adjusted_likelihood.idxmax()].adjusted_likelihood ])
        
        best_domain_properties_table_SET1.append(best_domain_properties)
    
    if write_output:
        best_domain_properties_table_SET1_DF = pd.DataFrame(best_domain_properties_table_SET1_table, columns=['gene', 'gene_size_kbp', 'extent_left_kbp', 'extent_right_kbp','L'])
        best_domain_properties_table_SET1_DF.to_csv( today+'table_of_domains'+'_set'+str(line_no)+'_chr'+the_chr+'.csv', sep='\t', index=False)
        
        gene_DHS_rejection_tableDF = pd.DataFrame(gene_DHS_rejection_table, columns=['gene', 'DHSind'])
        gene_DHS_rejection_tableDF.to_csv( today+'gene_DHS_rejections_set'+str(line_no)+'_chr'+the_chr+'.csv', index=False, sep='\t')
        
        gene_mixture_tableDF = pd.DataFrame(gene_mixture_table, columns=['gene'] + ['C'+str(i) for i in range(1,17)])
        gene_mixture_tableDF.to_csv( today+'gene_mixture_list'+str(line_no)+'_chr'+the_chr+'.csv', index=False, sep='\t')

        gene_performance_tableDF = pd.DataFrame(gene_performance_table, columns=['gene', 'W_size', 'W_coherence', 'W_foreground', 'W_centeredness'])
        gene_performance_tableDF.to_csv( today+'gene_performance_list'+str(line_no)+'_chr'+the_chr+'.csv', index=False, sep='\t')


        gene_DHS_ownership_tableDF = pd.DataFrame(gene_DHS_ownership_table, columns=['gene', 'DHSind'])
        gene_DHS_ownership_tableDF.to_csv( today+'gene_DHS_inds_set'+str(line_no)+'_chr'+the_chr+'.csv', index=False, sep='\t')
        

    allDHS_indices_SET1 = []
    for entry in best_domain_properties_table_SET1:
        if entry['chrom'] == 'chr'+chrom:
            allDHS_indices_SET1 += entry['final_DHS_indices_owned'].tolist()

    allcounts = pd.Series(allDHS_indices_SET1).value_counts()
    #uniqueness_ratio =  len(allcounts)/allcounts.sum() 
    #assignment_ratio = len(allcounts) / len(DHS_X)
    
    total_assigned_DHSs.append(allcounts.sum())
    unique_assigned_DHSs.append(len(allcounts))
    total_DHS_per_chrom.append(len(DHS_X))


uniqueness_ratio = np.sum(unique_assigned_DHSs) / np.sum(total_assigned_DHSs)
assignment_ratio = np.sum(unique_assigned_DHSs) / np.sum(total_DHS_per_chrom)

f = open(today+'_sweep_results_total_allchr.txt', 'a')
print(line_no,thef1, thef2, thef3, thef4, thef5, thecosthresh, uniqueness_ratio, assignment_ratio, file=f)
f.close()

