# Vocabulary

A framework for decomposing DNase I Hypersensitive Sites (DHSs) measured across 100s of biosamples using Non-Negative Matrix Factorization (NMF),
with the goal of developing a regulatory "Vocabulary" for the accessible human genome.

## Implementation

Object-oriented (OO) implementation, requiring Python 3 with `scikit-learn`, `numpy`, `scipy`, `matplotlib`, `pandas` installations.
The primary OO-developed code are distributed in the root directory.
Various analysis and visualization efforts are demonstrated in the `notebooks` directory using Jupyter notebooks.
Requisite data are stored in the `data` directory, when possible.
Several larger files that are incompatible with Github are available through [Zenodo](https://doi.org/10.5281/zenodo.3752359).
Supporting scripts to obtain peak or signal data for new unseen biosamples are available in `scripts`.

Authors: Alexander (Sasha) Muratov & Wouter Meuleman

## Notebooks
1. **Basic procedure using random binary data.**
Shows how to decompose a random binary matrix and measure reconstruction quality.
2. **Unsupervised metrics for choosing number of NMF components.**
Shows how we chose k=16 NMF components using unsupervised metrics.
3. **Demonstrating decomposition results of ENCODE DNase-Seq biosamples.**
Demonstrates how we decomposed a binary 733x3.6M matrix into a Vocabulary of regulatory components
4. **Visualizing decomposition results of ENCODE DNase-Seq biosamples.**
Visualizes the results of the decomposition in the previous notebook. Includes UMAP projection.
5. **Demonstrating the embedding of new biosamples into the coordinate space.**
Demonstrates how one can embed new (unseen) DNase-seq data into an existing NMF decomposition.
6. **Visualizing the embedding of new biosamples into the coordinate space.**
Visualizes the newly added biosamples alongside the existing 733. Includes UMAP projections.


