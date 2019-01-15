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


def get_barsortorder_OfficialOrder(relevantMatrix):
    # much more ugly but gets the job done
    # assumes rows are the data, columns are NMF components
    WSO = np.array([7,5,15,9,12,14,3,8,13,2,4,6,16,11,10,1]).astype(int) - 1
    WinningComponent = np.argmax(relevantMatrix, axis=1)
    barsortorder = np.array([])
    for i in range(relevantMatrix.shape[1]):
        desired_order = np.argsort(-relevantMatrix[:,WinningComponent[WSO[i]]])
        relevant_cut = WinningComponent[desired_order]==WSO[i]
        barsortorder = np.append(barsortorder, desired_order[relevant_cut])
    barsortorder = barsortorder.astype(int)
    return barsortorder
	
