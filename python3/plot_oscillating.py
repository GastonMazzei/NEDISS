import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

with open('/home/m4zz31/cppprojct/graphic/timeseries/result.pkl','rb') as f:
  data = pickle.load(f)


dataKeys = list(data.keys())
L = len(data[dataKeys[0]])
D = len(str(L))
MIN,MAX = 1,-1
for k in dataKeys: 
	assert(len(data[k])==L)
	if min(np.sin(data[k])) < MIN:
		MIN = min(np.sin(data[k]))
	if max(np.sin(data[k])) > MAX:
		MAX = max(np.sin(data[k]))

for i in range(L):
	v = []
	for k in dataKeys:
		v += [np.sin(data[k][i])]
	plt.bar(range(len(v)), v, color='r')
	plt.ylim(MIN-0.1,MAX+0.1)
	plt.savefig(f'graphic/timeseries/oscillation-files/{str(i).zfill(D)}.png')
	plt.clear()
	plt.close()


os.system(f'ffmpeg  -i "graphic/timeseries/oscillation-files/%0{D}d.png" -vcodec mpeg4 "graphic/video/result-oscillating.avi"')
