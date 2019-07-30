# EncodeNMF

Object-oriented implementation of the of the framework of methods I've developed to analyze results from decompositions of Non-negative matrix factorization (NMF).
The NMF itself is done with `scikit-learn` routines, but all the data manipulation, visualziation, interpretation requires extra code. 

the primary OO-developed code are distributed in the root directory, while all the scripts that utilize these routines are in `scripts`. Further analysis and data exploration for our component mixture project are found in the `notebooks` directory. Finally, old versions of certain routines are stored in `old` for archiving. 

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
- OONMF usage basics - decomposing a random binary matrix and measuring reconstruction quality .ipynb : this shows basic usage of my NMF package on random data. 
- OONMF analyzing results from ENCODE DNase-Seq .ipynb - this takes us through the process of visualizing the NMF decomposition with 16 components on Presence/Absence matrix for 733-sample masterlist. Also shows UMAP visualization and projection of new pancreatic cancer samples into the space.
- OONMF how we chose 16 components - unsupervised metrics.ipynb - this notebook shows some analysis of unsupervised metrics that we used in our decision making process for 16 components.
- OONMF using metadata for semi-supervised selection of k .ipynb - this notebook shows the additional constraints we obtain by trying to statistically associate NMF components with metadata labels for the biosamples
- OONMF the relationship between components, samples per DHS, and systems per DHS - this notebook takes a look at the DHS x Component matrix, and explores how the number of non-zero NMF component loadings per DHS correlates with other things.

# TODO 
- More technical stuff
  - PCA vs NMF
  - KL divergence vs L2 norm for NMF objective function
  - comparing data to random noise matrix
- comparing all the different groups and calculating statistics about them - DONE 
- setting up & analyzing GREAT runs 
- fisher tests for motif pairs in G2 components
- 2.5 D NMF - analysis of motif components
- Regulatory domains done 3 ways
  - "Regulatory Potential"
  - GMM / EM approach (TBD) 
  - Likelihood function based on chosen constraints
- comparison of reg domains to gene expression. 
