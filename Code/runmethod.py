from DeltaTrimax import DeltaTrimax
import numpy as np
import pandas as pd

#import data
data = pd.read_csv("../Data/data.xls",index_col=0)

#normalization
def logmean_centering(data):
	lx = data.apply(np.log10)
	a = np.mean(lx)
	b = lx-a
	N_data = b
	return N_data

#separate data base condition
k1 = N_data.iloc[:,0:10].values
k2 = N_data.iloc[:,10:20].values
k3 = N_data.iloc[:,20:30].values
k4 = N_data.iloc[:,30:40].values

#combine to 3d array
D = np.array([k1,k2,k3,k4])

#variavel
dt=float(sys.argv[1])
lm=float(Sys.argv[2])
#method
a = DeltaTrimax(D)
ga,ka,wa, msra, tqi = a.fit(dt,lm,n_triclusters=0)
#save output
np.savetxt("gen.txt",g,fmt="%0.f")
np.savetxt("kondisi.txt",k,fmt="%0.f")
np.savetxt("waktu.txt",w,fmt="%0.f")
np.savetxt("msr.txt",msr)
np.savetxt("tqi.txt",tqi)
