# Commands to run NEDISS


mpirun -oversubscribe -x SOLVER=1 -x TOPOLOGY=0 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct > tmp/trash.txt ; python3 trash_counter.py

