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
a.make_stacked_bar_plot(a.BasisD, a.Basis.T, today+'ATesto1.png', barsortorder=farsortorder)

randomdata = acquire_data()
PRC = a.precision_recall_curve(data=randomdata)
print(PRC.sort_values(by='F1')[-5:])

