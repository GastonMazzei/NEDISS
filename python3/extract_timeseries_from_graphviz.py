import os
import matplotlib.pyplot as plt
import pickle


DIR = 'graphic/program-output'
files = [x for x in os.listdir(DIR) if x[-3:]=='dot']


Ntot = 0
with open(f'{DIR}/{files[0]}','r') as f:
	for x in f.readlines():
		if '[label' in x:
			temp = int(x.split('n')[1].split('[')[0])
			if temp > Ntot:
				Ntot = temp;

v = {i:[0 for _ in range(len(files))] for i in range(Ntot+1)}
for fn in os.listdir(DIR):
	with open(f'{DIR}/{fn}','r') as f:
		i = int(fn.split('.')[1])
		for x in f.readlines():
			if '[label' in x:
				v[int(x.split('n')[1].split('[')[0])][i] = float(x.split('=')[1].split(']')[0])

		

for k in v.keys():
	plt.plot(v[k])

plt.savefig('graphic/image/result.png')

with open('graphic/timeseries/result.pkl','wb') as f:
	pickle.dump(v,f)
