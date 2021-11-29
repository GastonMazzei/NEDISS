import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define two PI
TWOPI = np.pi * 2


# Open the simulation configuration file and get relevant info
with open('simulation.conf.sh','r') as f:
	for x in f.readlines():
		if 'h=' in x and x[0]=='h':
			h = float(x.split('=')[1][:-2])
		if 'J=' in x and x[0]=='J':
			J = float(x.split('=')[1][:-2])
		if 'NRUNS=' in x:
			NRUNS = float(x.split('=')[1][:-2])
		if 'EQNUMBER=' in x:
			EQNUMBER = float(x.split('=')[1][:-2])
		if 'SAMPLING_FREQ=' in x:
			SAMPLING_FREQ = float(x.split('=')[1][:-2])
		if 'WMIN=' in x:
			WMIN = float(x.split('=')[1][:-2])
		if 'WMAX=' in x:
			WMAX = float(x.split('=')[1][:-2])
		if 'TOPOLOGY=' in x:
			TOPOLOGY = float(x.split('=')[1][:-2])
		if 'proba=' in x:
			proba = float(x.split('=')[1][:-2])
		if 'kneigh=' in x:
			kneigh = float(x.split('=')[1][:-2])

# MEANWHILE, AS 'h' IS NOT CURRENTLY PASSED AS PARAMETER WE PUT IT MANUALLY
h = 0.01

def func(x, T, a, b):
    return a * np.exp(-x / T) + b


with open('graphic/timeseries/result.pkl','rb') as f:
	data = pickle.load(f)

L = max(data.keys())+1

v = np.asarray([data[k_] for k_ in range(L)])


def compute_distance(a,b):
	return min( [  abs((a-b) % TWOPI) , abs((b-a) % TWOPI) ] )	

def get_max_distance(y):
	initial_d = -np.inf
	Len = len(y)
	for i in range(Len):
		for j in range(i+1, Len-i-1):
			local_d = compute_distance(y[i],y[j])
			if local_d > initial_d:
				initial_d = local_d	
	return initial_d

dists = np.asarray([get_max_distance(v[:,i]) for i in range(v.shape[1])])

ETOL = abs(max(dists) - min(dists))/v.shape[1] 
mask = np.where (np.abs(np.diff(dists))>ETOL,True,False)
new_dists = dists[1:][mask]

print(f'sampling freq is {SAMPLING_FREQ}')
xdata = np.asarray([h * x for x in range(len(new_dists))])
ydata = new_dists
popt, pcov = curve_fit(func, xdata, ydata, 
			#bounds=(0, [3., 1., 0.5]) # currently no bounds.
				)

plt.plot(xdata, func(xdata, *popt), 'g--', lw=3, alpha=0.8, label='fit: T=%9.4f, a=%4.2f, b=%4.2f' % tuple(popt))
plt.scatter(xdata, new_dists, lw=3, c='b', ls=':', label='data')
plt.legend()
plt.title(r'Decay towards a locked state fitted with $ae^{\frac{-t}{T}}+b$')
plt.savefig('graphic/image/maxdistovertime.png')


T = popt[0] * SAMPLING_FREQ
#	    T, h, J, NRUNS, EQNUMBER, SAMPLING_FREQ, WMIN, WMAX, TOPOLOGY, proba, kneigh
myCsvRow = f'{T},{h},{J},{NRUNS},{EQNUMBER},{SAMPLING_FREQ},{WMIN},{WMAX},{TOPOLOGY},{proba},{kneigh}\n'
with open('graphic/timeseries/total-results.csv','a') as fd:
    fd.write(myCsvRow)


