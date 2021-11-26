import numpy as np
import os

filenames = os.listdir('graphic/program-output')
TOT = len(filenames)

d = [[] for _ in range(TOT)]
for i in range(TOT):
	with open(f'graphic/program-output/graphviz.{i}.dot','r') as f:
		d[i] = f.readlines()


# Compute a range 
vals = []
for i in range(TOT):
	for l in d[i]:
		if '[label' in l:
			vals += [float(l.split('[')[1].split('=')[1].split(']')[0])]
MIN, MAX = min(np.sin(vals)), max(np.sin(vals))


# Global values
COLORS = ['#0d0887', '#46039f', '#7201a8', '#9c179e', '#bd3786', '#d8576b', '#ed7953', '#fb9f3a', '#fdca26', '#f0f921']
L = len(COLORS)
VALUES = np.linspace(MIN,MAX,L)


def ftoh(v):
	"""
	convert value to hex
	"""
	return COLORS[np.argmin(np.abs(VALUES - np.sin(v)))] 


for i in range(TOT):
	with open(f'graphic/processed-program-output/graphviz.{i}.dot','w') as f:
		for l in d[i]:
			if '->' in l:
				f.write(l[:-2] + ' [dir=none];\n')
			elif '[label' in l:
				[first,second] = l.split('[')
				v = float(second.split('=')[1].split(']')[0])
				f.write(first+f'[style=filled, fillcolor="{ftoh(v)}"];\n')
			else:
				f.write(l)

