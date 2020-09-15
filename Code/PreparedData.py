import pandas as pd
import numpy as np

class PreparedData():
	def data():
		data = pd.read_csv("Data/data.xls",index_col=0)
		return data

	#Normalisasi Logmean-Centerin
	def logmean_centering(data):
		lx = data.apply(np.log10)
		a = np.mean(lx)
		b = lx-a
		N_data = b
		return N_data
	def fit():
		data = data()
		N_data = logmean_centering(data)
		#separate data base condition
		k1 = N_data.iloc[:,0:10].values
		k2 = N_data.iloc[:,10:20].values
		k3 = N_data.iloc[:,20:30].values
		k4 = N_data.iloc[:,30:40].values

		return k1, k2, k3, k4


#compound data as 3d array condition-gene-time
