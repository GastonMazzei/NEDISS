import os, sys
import numpy as np


def command_maker(S):
	THR = np.random.randint(3,10)
	P = np.random.randint(3,10)
	N = np.random.randint(P*THR*2, P*THR*10)
	command = f"mpirun -x NNODES={N} -x TEST=1 -x SEED={S} -x OMP_THREAD_LIMIT={THR}  -x  OMP_NESTED=true -n {P}  cmake-build-debug/cppprojct"
	return THR,P,N,command

errors = 0
success = 0

while (errors == 0):
	THR,P,N,command = command_maker(np.random.randint(0,999999))
	response = os.system(command + " >  tmp/trash.txt")
	if response == 0: success += 1
	else: errors += 1
	try: os.remove('tmp/trash.txt')
	except: pass
	if (response != 0): print(f'This script has found its death :-0. This case has Nthr:{THR} Nproc:{P} Nnodes:{N}')
	else: print(f'So far no errors after {success} runs ;-). This case has Nthr:{THR} Nproc:{P} Nnodes:{N}')

print(f'\n--------------\nFinal report: observed {errors} errors in {success+errors} runs\n------------------\n')
