ClusterMode = True

import datetime
import time
from datetime import date
import numpy as np
today = str(date.today())
if (ClusterMode):
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

def mytime():
    ts = time.time()
    sttime = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H:%M:%S')
    return sttime

def increase_axis_fontsize(axis_fontsize=30):
    ax = plt.gca()
    ticklabels = ax.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(axis_fontsize)
        label.set_family('serif')
    ticklabels = ax.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(axis_fontsize)
        label.set_family('serif')

def get_today():
	today = str(date.today())
	return today

def get_barsortorder(relevantMatrix):
    # assumes rows are the data, columns are NMF components
    WinningComponent = np.argmax(relevantMatrix, axis=1)
    barsortorder = np.array([])
    for i in range(relevantMatrix.shape[1]):
        desired_order = np.argsort(-relevantMatrix[:,i])
        relevant_cut = WinningComponent[desired_order]==i
        barsortorder = np.append(barsortorder, desired_order[relevant_cut])
    barsortorder = barsortorder.astype(int)
    return barsortorder
    Basis_names = pd.read_table('colTitles.651samples_friendly.txt', header=None)[1].values


	