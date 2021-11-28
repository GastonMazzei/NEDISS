# Commands to run NEDISS


mpirun -oversubscribe  -x NRUNS=2 -x EQNUMBER=0 -x SOLVER=1 -x SOLVERDEPTH=2 -x K1=0.1666666666666666 -x K2=0.33333333333333333 -x K3=0.333333333333333 -x K4=0.16666666666666  -x TOPOLOGY=0 -x kneigh=3 -x proba=0.1234  -x NNODES=12 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=5  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct


