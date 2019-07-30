import numpy as np
import pandas as pd



def acquire_data():
    # fill in for whatever application
	pass
	#print('starting read at ', mytime(), flush=True)
	A = pd.read_table('masterListPresenceAbsenceMatrix.651samples.FDR0.0010.hg38.bed', header=None)
	#print('finished read at ', mytime(), flush=True)
	return A.drop([A.columns[0], A.columns[1], A.columns[2]], axis=1).values.T


C = acquire_data()

f = open('sparse651.txt', 'w')

print(C.shape[0], C.shape[1], file=f)

for i in range(C.shape[0]):
  for j in range(C.shape[1]):
   if C[i][j] > 0:
    print(i,j,C[i][j],file=f)
    
f.close()
