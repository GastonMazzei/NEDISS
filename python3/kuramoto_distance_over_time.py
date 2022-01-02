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
			h = float(x.split('=')[1].split(';')[0])
		if 'J=' in x and x[0]=='J':
			J = float(x.split('=')[1].split(';')[0])
		if 'NRUNS=' in x:
			NRUNS = float(x.split('=')[1].split(';')[0])
		if 'EQNUMBER=' in x:
			EQNUMBER = float(x.split('=')[1].split(';')[0])
		if 'SAMPLING_FREQ=' in x:
			SAMPLING_FREQ = float(x.split('=')[1].split(';')[0])
		if 'WMIN=' in x:
			WMIN = float(x.split('=')[1].split(';')[0])
		if 'WMAX=' in x:
			WMAX = float(x.split('=')[1].split(';')[0])
		if 'TOPOLOGY=' in x:
			TOPOLOGY = float(x.split('=')[1].split(';')[0])
		if 'proba=' in x:
			proba = float(x.split('=')[1].split(';')[0])
		if 'kneigh=' in x:
			kneigh = float(x.split('=')[1].split(';')[0])

# MEANWHILE, AS 'h' IS NOT CURRENTLY PASSED AS PARAMETER WE PUT IT MANUALLY
h = 0.01

def func(x, T, a, b):
    return a * np.exp(-x / T) + b


with open('graphic/timeseries/result.pkl','rb') as f:
	data = pickle.load(f)

L = max(data.keys())+1

v = np.asarray([data[k_] for k_ in range(L)])


def compute_distance(a,b):
	"""
	S1 distance in period 2*PI :-)
	"""
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

ixs = np.asarray(range(len(dists)-1))[mask]


# Flagged for kill after next debug 01 jan 2022
# Filter ixs to keep the last fully connected blob :-)
# Don't do it... disconnected is not bad :-)
#
#new_ixs = [ixs[-1]]
#for ix in ixs.tolist()[:-1][::-1]:
#	if ix == new_ixs[-1] - 1:
#		new_ixs.append(ix)
#	else:
#		break
#new_ixs.reverse()

new_ixs = ixs.copy()

# Compute the asymptotic value
Lnd = len(new_ixs)
L = len(dists)
Llast = max([1,int(L*0.1)])
asymptote = np.mean(dists[-Llast:])

# Go through a case-by-case analysis to select the exponential part
if np.mean(dists[new_ixs[:2]]) > np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]) and np.mean(dists[new_ixs[-2:]]) > np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]):
	# it's a  
	#   \    /
	#    \  /  
	#     \/  
	if abs(np.mean(dists[new_ixs[-2:]])-asymptote)<0.1*ETOL*v.shape[1]*2:
		# The full picture is probably
		#         _______ 
		#   \    /
		#    \  /  
		#     \/  
		# Thus we compute the minimum, i.e.
		minval = min(dists[new_ixs])
		ixmin = new_ixs[dists[new_ixs] == minval][0]
		# And we keep only that section
		new_dists = dists[1:][new_ixs[ixmin:].tolist()]# + list(range(new_ixs[-1]+1,len(dists)-1))]
	else:		
		# The full picture is not clear
		plt.plot(new_ixs, dists[new_ixs], c='r', label='selected')
		plt.plot(dists, c='b', label='full')
		plt.legend()
		plt.title('We had a problem detecting the exponential convergence section... can you please input where it starts manually? (aprox x values)')
		plt.show()		
		xs = input('Input the x value in integer format, e.g. "15"... Your answer: ')
		ixmin = int(xs)
		new_dists = dists[ixmin:]
elif np.mean(dists[new_ixs[:2]]) <= np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]) and np.mean(dists[new_ixs[-2:]]) <= np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]):
	# it's a  
	#    /\   
	#   /  \    
	#  /    \ 
	if abs(np.mean(dists[new_ixs[-2:]])-asymptote)<0.1*ETOL*v.shape[1]*2:
		# The full picture is probably
		#    /\   
		#   /  \    
		#  /    \_____ 
		# Thus we compute the maximum, i.e.
		maxval = max(dists[new_ixs])
		ixmax = new_ixs[dists[new_ixs] == maxval][0]
		# And we keep only that section
		new_dists = dists[1:][new_ixs[ixmax:].tolist()] # + list(range(new_ixs[-1]+1,len(dists)-1))]
	else:		
		# The full picture is not clear
		plt.plot(new_ixs, dists[new_ixs], c='r', label='selected')
		plt.plot(dists, c='b', label='full')
		plt.legend()
		plt.title('We had a problem detecting the exponential convergence section... can you please input where it starts manually? (aprox x values)')
		plt.show()		
		xs = input('Input the x value in integer format, e.g. "15"... Your answer: ')
		ixmin = int(xs)
		new_dists = dists[ixmin:]
elif np.mean(dists[new_ixs[:2]]) > np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]) and np.mean(dists[new_ixs[-2:]]) <= np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]):
	# it's a  
	#   \    
	#    \    
	#     \/\/\__  
	if abs(np.mean(dists[new_ixs[-2:]])-asymptote)<0.1*ETOL*v.shape[1]*2:
		# The full picture is probably
		#   \    
		#    \    
		#     \/\/\_________  
		# Thus we keep it all.
		new_dists = dists[1:][new_ixs.tolist()]# + list(range(new_ixs[-1]+1,len(dists)-1))]
	else:		
		# The full picture is not clear
		plt.plot(new_ixs, dists[new_ixs], c='r', label='selected')
		plt.plot(dists, c='b', label='full')
		plt.legend()
		plt.title('We had a problem detecting the exponential convergence section... can you please input where it starts manually? (aprox x values)')
		plt.show()		
		xs = input('Input the x value in integer format, e.g. "15"... Your answer: ')
		ixmin = int(xs)
		new_dists = dists[ixmin:]
elif np.mean(dists[new_ixs[:2]]) <= np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]) and np.mean(dists[new_ixs[-2:]]) >= np.mean(dists[new_ixs[Lnd // 2-1:Lnd // 2+1]]):
	# it's a ___ 
	#     /\/   
	#    /    
	#   /   
	if abs(np.mean(dists[new_ixs[-2:]])-asymptote)<0.1*ETOL:
		# The full picture is probably
		#        ___________
		#     /\/   
		#    /    
		#   /   
		# Thus we keep it all.
		new_dists = dists[1:][new_ixs.tolist()]# + list(range(new_ixs[-1]+1,len(dists)-1))]
	else:		
		# The full picture is not clear
		plt.plot(new_ixs, dists[new_ixs], c='r', label='selected')
		plt.plot(dists, c='b', label='full')
		plt.legend()
		plt.title('We had a problem detecting the exponential convergence section... can you please input where it starts manually? (aprox x values)')
		plt.show()		
		xs = input('Input the x value in integer format, e.g. "15"... Your answer: ')
		ixmin = int(xs)
		new_dists = dists[ixmin:]
else:
	print("Last")
	# The full picture is not clear
	plt.plot(new_ixs, dists[new_ixs], c='r', label='selected')
	plt.plot(dists, c='b', label='full')
	plt.legend()
	plt.title('We had a problem detecting the exponential convergence section... can you please input where it starts manually? (aprox x values)')
	plt.show()		
	xs = input('Input the x value in integer format, e.g. "15"... Your answer: ')
	ixmin = int(xs)
	new_dists = dists[ixmin:]

# Flagged for kill after next debug, 01 jan 2022
#
#ixs = np.asarray(range(len(dists)-1))[mask]
#new_dists = dists[1:][ixs.tolist() + list(range(ixs[-1]+1,len(dists)-1))]
#plt.plot(new_dists);plt.show()
#new_dists = dists[1:][mask]


# Function fitting :-)
print(f'[WARNING] Using sampling freq {SAMPLING_FREQ} (comes from config) and DT {h} (still not part of the config, it\'s hardcoded both here and in the compiled code)')
xdata = np.asarray([h * x for x in range(len(new_dists))])
ydata = new_dists
popt, pcov = curve_fit(func, xdata, ydata, 
			bounds=([0,-np.inf,0], [np.inf, np.inf, 3.14]) # currently no bounds.
				)
								# Short version: only T				# Long version: all coefficients :-)
plt.plot(xdata, func(xdata, *popt), 'g--', lw=3, alpha=0.8, label='fit: T=%9.4f' % popt[0]) #, label='fit: T=%9.4f, a=%4.2f, b=%4.2f' % tuple(popt))
plt.scatter(xdata, new_dists, lw=3, c='b', ls=':', label='data')
plt.legend()
plt.title(r'Decay towards a locked state fitted with $exp(-t/T)$'+'\n'+r'($\it{Max}$ $\it{Phase}$ $\it{difference}$ $\it{over}$ $\it{time,}$ $\it{The}$ $\it{limit}$ $\it{value}$ $\it{is}$: '+f'{round(popt[2],2)}')
plt.savefig('graphic/image/maxdistovertime.png')
# Output size: 640 x 480
# todo: put the entire graph and add a section of zoomed axis: https://matplotlib.org/stable/gallery/subplots_axes_and_figures/zoom_inset_axes.html
# actually we can double the vertical dimension and just vertically stack them :-)

T = popt[0] * SAMPLING_FREQ
#	    T, h, J, NRUNS, EQNUMBER, SAMPLING_FREQ, WMIN, WMAX, TOPOLOGY, proba, kneigh
myCsvRow = f'{T},{h},{J},{NRUNS},{EQNUMBER},{SAMPLING_FREQ},{WMIN},{WMAX},{TOPOLOGY},{proba},{kneigh}\n'
with open('graphic/timeseries/total-results.csv','a') as fd:
    fd.write(myCsvRow)


