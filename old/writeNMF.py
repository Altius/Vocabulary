import sys
import datetime
import time


writefiles = True


myrandom = 0

def mytime():
	ts = time.time()
	sttime = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H:%M:%S')
	return sttime
print('starting at ', mytime(), flush=True)


from sklearn.decomposition import NMF
import numpy as np
import pandas as pd
print('starting read at ', mytime(), flush=True)

A = pd.read_table('masterListPresenceAbsence_804samples_FDR0.0010_hg38.bed')
print('finished read at ', mytime(), flush=True)

a = A.drop([A.columns[0], A.columns[1], A.columns[2]], axis=1).values

a = a.T
print('finished drop at ', mytime(), flush=True)


answer = 8

model = NMF(n_components=answer, init='random', random_state=myrandom)
print('starting NMF at ', mytime(), flush=True)

newdaslog = model.fit_transform(a) 
print('done with NMF at ', mytime(), flush=True)
newdaslog2 = model.components_

if (writefiles):
	outfile = today+'NMF_Ncomps'+str(answer)+'daslog.npy'
	np.save(outfile, newdaslog)	
	outfile = today+'NMF_Ncomps'+str(answer)+'daslog2.npy'
	np.save(outfile, newdaslog2)	

