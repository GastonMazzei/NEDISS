# Commands to run NEDISS

# TEST=2 requires NRUNS=50 (or other number different than 50!)
mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=0 -x SOLVERDEPTH=3 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=0 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct

