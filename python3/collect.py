with open('results.txt','r') as f: d=f.readlines()
L = [int(d[i]) for i in range(len(d)) if i%3==0]
tseq = [float(d[i]) for i in range(len(d)) if i%3==1]
tpar = [float(d[i]) for i in range(len(d)) if i%3==2]

import matplotlib.pyplot as plt

plt.plot(L,tseq, label='Tseq');plt.plot(L,tpar,label='Tpar');plt.legend();plt.show()
