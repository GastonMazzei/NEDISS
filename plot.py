import matplotlib.pyplot as plt
import os
import numpy as np

# Define processor
def process(d: str) -> list:
    """
    Recieves a string and produces a list of floats
    with all the comma-delimited numerical values
    that were in the string
    """
    y = []
    for x in d.split(','):
        try:
            y += [float(x)]
        except:
            print(f'[info] this value could not be casted as float: {x}')
    return y

# Collect all the oscillator's result names
nms = [x for x in os.listdir('cmake-build-debug') if 'osc' in x]

# Define a dictionary and load the data
d = {}
for x in nms:
    with open(f'cmake-build-debug/{x}', 'r') as f:
        temp = f.readlines()[0]
    ix = x.split('osc')[1].split('.')[0]
    d[ix] = process(temp)


# Plot the data for each oscillator
print(d.keys())
if [False, True][1]:
    for nm,v in d.items():
        plt.plot(v, label=nm)
        #plt.plot(np.sin(v), label=nm)
else:
    plt.plot(np.asarray(d['1'])-np.asarray(d['2']))
plt.legend()
plt.show()
