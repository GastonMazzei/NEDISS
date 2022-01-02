
import os

try: os.remove('graphic/timeseries/total-results.csv')
except: pass

# Define a configuration updater :-)
def update_config(N,P,K,J,T):
    with open('simulation.conf.sh','r') as f:
        d = f.readlines()
    for i,x in enumerate(d):
        if 'TOPOLOGY' in x:
            d[i] = f'TOPOLOGY={T};\n'
        if 'NNODES' in x:
            d[i] = f'NNODES={N};\n'
        if 'proba' in x:
            d[i] = f'proba={P};\n'
        if 'kneigh' in x:
            d[i] = f'kneigh={K};\n'
        if 'J=' in x or 'J =' in x:
            d[i] = f'J={J};\n'
    with open('simulation.conf.sh','w') as f:
        for x in d:
            f.write(x)

# User can set the range here
try:
    from period_analysis_configuration import SAMPLES,TOPOLOGY,NNODES,PROBA,KNEIGH,J
    print("Successfully loaded the period_analysis_configuration file!")
except:
    print("Failed to load the period_analysis_configuration file!")
    SAMPLES = 1
    TOPOLOGY = 0
    NNODES =  [20 + i * 10 for i in range(20)]
    PROBA = [0.5]
    KNEIGH = [3]
    J = [2]

# Compute lengths and assert that only one will ve varied
LNNODES = len(NNODES)
LPROBA = len(PROBA)
LKNEIGH = len(KNEIGH)
LJ = len(J)
dfz = False
Lmax = 1
for x in [LNNODES, LPROBA, LKNEIGH, LJ]:
    if x!=1:
        if x>Lmax: Lmax = x
        if dfz: raise Exception("[ERROR] Only one variable can be varied during the exploration")
        else: dfz = True
print("Successfully passed the filtering criteria :-)")


counter = 0
Lmax *= SAMPLES
for _ in range(SAMPLES):
    for N in NNODES:
        for P in PROBA:
            for K in KNEIGH:
                for _J in J:
                    update_config(N,P,K,_J,TOPOLOGY)
                    print(f'\n*******\nStarting {counter+1} of {Lmax}\n*******\n')
                    os.system('/bin/bash ./capture_period.sh > temp.txt')
                    os.remove('temp.txt')
                    os.system('python3 python3/store_periodsweep_data.py')
                    counter += 1
