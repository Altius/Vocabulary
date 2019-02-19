# EncodeNMF

Object-oriented implementation of the of the framework of methods I've developed to analyze results from decompositions of Non-negative matrix factorization (NMF).
The NMF itself is done with `scikit-learn` routines, but all the data manipulation, visualziation, interpretation requires extra code. 

the primary OO-developed code are distributed in the root directory, while all the scripts that utilize these routines are in `usage`. 
Because I am not always good, I also have a non-object-oriented implementation of some routines in `imparative`. A bunch of data exploration has been done in jupyter notebooks in the `notebooks` directory. Finally, old versions of certain routines are stored in `old` for archiving. 

To run the code as we did in the ENCODE 3 masterlist paper witk `k=16`, simply download the repository and run:

```
python OONMF_compute_presence_NNDSVD_O.py
```

to run multiple iterations at once, instead run 

```
python OONMF_multirun.py
```

Note: requires Python 3 with `scikit-learn`, `numpy`, `scipy`, `matplotlib`, `pandas` installations.
