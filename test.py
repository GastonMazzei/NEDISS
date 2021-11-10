import time
import os,sys
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(0) # SET THE SEED!


def runscript(timeout, thread_limit, number_processes):
	SEED = int(np.random.randint(0,99999999))
	return f"mpiexec --timeout {int(timeout)} -x SEED={SEED} -x OMP_THREAD_LIMIT={thread_limit}  -x  OMP_NESTED=true -n {number_processes}  cmake-build-debug/cppprojct > temp/trash.txt"


if __name__=='__main__':
	try: os.remove('temp/trash.txt')
	except: pass

	NPROCS, NTHREADSMAX, NRUNS = 8, 10, 20
	NPROCS_v, NTHREADSMAX_v = [],[]

	# measure characteristic time 
	CHARACTERISTIC_TIME = 120 #30 SECONDS initial guess
	#start = time.time()
	output = os.system(runscript(CHARACTERISTIC_TIME, NTHREADSMAX, NPROCS))
	end = time.time()
	#CHARACTERISTIC_TIME = end - start + 5
	if (output == 0): print(f'output was: {output}')
	else: 
		print('there was an error')
		sys.exit(output)

	# Now that we have the characteristic time,
	# run the script several times to test for memory errors.
	times = []
	for i in range(NRUNS):
		NPROCS = np.random.randint(2,8)
		NPROCS_v += [NPROCS].copy()
		NTHREADSMAX = np.random.randint(3,10)
		NTHREADSMAX_v += [NTHREADSMAX].copy()
		start = time.time()
		output = os.system(runscript(CHARACTERISTIC_TIME, NTHREADSMAX, NPROCS))

		if (output == 0): 
			print(f'output was: {output}')
			end = time.time()
			times.append(end-start)

		elif (output == 256): 
			print (f'we ran out of time, but it wasnt sigfault. lap {i+1} of {NRUNS}, output was: {output} '+
						f'where Nprocs and Nthreadsmax were {NPROCS} & {NTHREADSMAX} respectively.')
			*NTHREADSMAX_v,_ = NTHREADSMAX_v
			*NPROCS_v,_ = NPROCS_v
		else: 
			print(f'there was an error at lap {i+1} of {NRUNS}, output was: {output}')
			print(f'Nprocs and Nthreadsmax were {NPROCS} & {NTHREADSMAX} respectively.')
			print(f'I should have died at lap {i+1}.')
			print('WARNING: POTENTIAL SIGFAULT FOUND')
			*NTHREADSMAX_v,_ = NTHREADSMAX_v
			*NPROCS_v,_ = NPROCS_v
			#sys.exit(output)
	
	# Final report:
	try: os.remove('temp/trash.txt')
	except: pass

	print(f'There was no exit status, and in average the script took {np.mean(times)} seconds for Nprocs: {NPROCS} and NthreadsMax: {NTHREADSMAX}')

	f,ax = plt.subplots(1,3,figsize=(15,5))
	ax[0].scatter(*zip(*sorted([(NTHREADSMAX_v[i], times[i]) for i in range(len(times))], key=lambda tup: tup[1])), c='r')
	ax[1].scatter(*zip(*sorted([(NPROCS_v[i], times[i]) for i in range(len(times))], key=lambda tup: tup[1])), c='b')
	ax[0].plot(*zip(*sorted([(NTHREADSMAX_v[i], times[i]) for i in range(len(times))], key=lambda tup: tup[1])), c='r')
	ax[1].plot(*zip(*sorted([(NPROCS_v[i], times[i]) for i in range(len(times))], key=lambda tup: tup[1])), c='b')
	ax[0].set_title('Time vs Threads (for potentially different Nprocs in each)')
	ax[1].set_title('Time vs Procs (for potentially different Nthreads in each)')
	ax[2].scatter(*zip(*sorted([(NPROCS_v[i]*NTHREADSMAX_v[i], times[i]) for i in range(len(times))], key=lambda tup: tup[1])), c='k')
	ax[2].plot(*zip(*sorted([(NPROCS_v[i]*NTHREADSMAX_v[i], times[i]) for i in range(len(times))], key=lambda tup: tup[1])), c='k')
	ax[2].set_title('Time vs Threads*Procs')
	plt.show()

