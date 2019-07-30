import numpy as np
import pandas as pd
import OONMF
from OONMFhelpers import *

today  = get_today()
def acquire_data():
    # fill in for whatever application
    pass
    
    data = np.random.rand(200, 200)
    data = (data>0.5).astype(int)
    return data
    
a = OONMF.NMFobject()
myMatrix = np.random.rand(200, 200)
a.performNMF(myMatrix)
farsortorder = get_barsortorder(a.Basis)
#a.make_stacked_bar_plot(a.BasisD, a.Basis.T, today+'ATesto1.png', barsortorder=farsortorder)

randomdata = acquire_data()

a.NormedMixture =   a.Mixture / np.sum(a.Mixture, axis=0)
a.NormedBasis =   (a.Basis.T / np.sum(a.Basis.T, axis=0)).T



a.make_stacked_bar_plot(a.BasisD, a.NormedBasis.T, today+'normedBasisRandom'+str(a.Ncomps)+'.png', barsortorder=farsortorder, names=a.Basis_Names)

a.make_stacked_bar_plot(a.MixtureD, a.NormedMixture, today+'normedMixtureRandom'+str(a.Ncomps)+'.png')

