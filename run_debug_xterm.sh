# Commands to run NEDISS



mpirun -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -n 3 xterm -e  gdb  -ex run  cmake-build-debug/cppprojct
