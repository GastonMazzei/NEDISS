# Commands to run NEDISS

mpirun -x NNODES=20 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=5  -x  OMP_NESTED=true -n 5  cmake-build-debug/cppprojct

