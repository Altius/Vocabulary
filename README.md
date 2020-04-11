# EncodeNMF

Object-oriented implementation of the framework of methods developed to analyze results from decompositions of Non-negative matrix factorization (NMF).
The NMF itself is done with `scikit-learn` routines, but all the data manipulation, visualization, interpretation requires extra code. 

the primary OO-developed code are distributed in the root directory, while all the scripts that utilize these routines are in `scripts`. 
Further analysis and data exploration for our component mixture project are found in the `notebooks` directory. Finally, old versions of certain routines are stored in `old` for archiving. 

To run the code as we did in the ENCODE 3 masterlist paper witk `k=16`, simply download the repository and run:

```
python scripts/OONMF_compute_presence_NNDSVD_O.py
```

to run multiple iterations at once, instead run 

```
python scripts/OONMF_multirun.py
```

Note: requires Python 3 with `scikit-learn`, `numpy`, `scipy`, `matplotlib`, `pandas` installations.

# Notebooks of interest:
- Chapter 1 usage basics - decomposing a random binary matrix and measuring reconstruction quality .ipynb : this shows basic usage of my NMF package on random data. 
- Chapter 2 analyzing results from ENCODE DNase-Seq .ipynb - this takes us through the process of visualizing the NMF decomposition with 16 components on Presence/Absence matrix for 733-sample masterlist. Also shows UMAP visualization and projection of new pancreatic cancer samples into the space.
- Chapter 3 how we chose 16 components - unsupervised metrics.ipynb - this notebook shows some analysis of unsupervised metrics that we used in our decision making process for 16 components.
- Chapter 4 using metadata for semi-supervised selection of k .ipynb - this notebook shows the additional constraints we obtain by trying to statistically associate NMF components with metadata labels for the biosamples
- Chapter 5 the relationship between components, samples per DHS, and systems per DHS - this notebook takes a look at the DHS x Component matrix, and explores how the number of non-zero NMF component loadings per DHS correlates with other things.
- Chapter 6 dividing DHSs into groups based on their component mixture.ipynb - this notebook introduces the idea of dividing the DHSs into "chromatic groups". It shows which component combinations are most common for DHSs that consist of 1, 2, or more different NMF components.
- Chapter 7 visualize motif cluster projection.ipynb - Using a method we call "2.5D NMF", we project regulatory factor motifs into the coordinate system defined by our biosample x DHS matrix. Here, we visualize these motifs in our coordinate system, and examine the reconstruction accuracy. We also compare with motif enrichment.
- Chapter 8 understanding the relationship between projected motifs and DHSs.ipynb - here we explore the possibility of using the projected motif component loadings to find the DHSs that actually contain those motifs. It works ok but not great.
- Chapter 9 This notebook demonstrates how I set up GREAT runs 

Here are a few technical notebooks i haven't had as much time to annotate
these can be found in  `notebooks/technical/` 



## Regulatory Domains:
- 8-15-19 check best model for all chr.ipynb this shows how we selected optimal paramters for the maximum likelihood model
- 8-15-19 new parameter grid search.ipynb this notebook documents how i used the maximum likeliihood model to find the boundaries of the regulatory domain of genes
- 9-03-19 fiducial results of parameter search.ipynb this notebook shows how i created the "Gene DHS masterlist" file

## gene-pathway connections (GREAT)

- 5-22-19 PCA on gene pathway pairs.ipynb  PCA analysis on gene-pathway pairs 
- 5-31-19 comprehensive pathway comparison .ipynb which pathways are really specific ot G2 DHSs?
- 6-06-19 combine enriched processes.ipynb "venn diagram" style list showing enriched pathways exclusive to and shared between most important component pairs
- 6-18-19 run the GREAT they hate G2.ipynb this is how I run GREAT using the R package rGREAT to check for enrichment


## TF motif analysis (not 2.5 D NMF version)

- 1-23-19 fisher exact tests.ipynb fisher exact test to show statistical significance of motif pairs co-occuring in DHSs of specific NMF component combinations
- 2-06-19 attempt to work with fimo - validation 7-30-19.ipynb  notebook that shows how i made a matrix for Presence/Absence of FIMO motifs used in other analysis
- 2-14-19 clean FIMO-based representation.ipynb - shows FIMO motif analysis figures for component / motif connection
- 3-04-19 motif vs motif significance simplified.ipynb which motif pairs are enriched in G2? From among a list of the top ~80 most enriched motifs overall


## Miscellaneous 
- 2-20-19 mean signal and N samples for Gs.ipynb shows mean signal in DHSs for G1, G2, G3, G4 component pair combinations
- 7-24-19 NMF on new CD3+ dataset.ipynb full analysis only on CD3+ stuff
- OONMF more metadata analysis.ipynb extra metadata categories screened







# TODO  
- More technical stuff
  - PCA vs NMF - can produce upon request
  - KL divergence vs L2 norm for NMF objective function - can produce upon request
  - comparing data to random noise matrix - can produce upon request
- comparing all the different groups and calculating statistics about them - DONE 
- setting up & analyzing GREAT runs  - DONE in technical notebooks
- fisher tests for motif pairs in G2 components - DONE in technical notebooks
- 2.5 D NMF - analysis of motif components - DONE
- Regulatory domains done 3 ways
  - "Regulatory Potential" - DONE in technical notebooks ( I believe)
  - GMM / EM approach (TBD)  
  - Likelihood function based on chosen constraints (prioritizing this) (DONE)
- comparison of reg domains to gene expression. (needs to be done, probably not by me)
