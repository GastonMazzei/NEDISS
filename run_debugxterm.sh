# Commands to run NEDISS




source simulation.conf.sh

mpirun -oversubscribe  -x SAMPLING_FREQ=$SAMPLING_FREQ -x EQNUMBER=$EQNUMBER -x SOLVER=$SOLVER -x SOLVERDEPTH=$SOLVERDEPTH -x K1=$K1 -x K2=$K2 -x K3=$K3 -x K4=$K4  -x TOPOLOGY=$TOPOLOGY -x J=$J -x WMIN=$WMIN -x WMAX=$WMAX -x kneigh=$kneigh -x proba=$proba -x NRUNS=$NRUNS -x NNODES=$NNODES -x TEST=1 -x SEED=$SEED -x OMP_THREAD_LIMIT=$OMP_THREAD_LIMIT  -x  OMP_NESTED=true -np $NP  xterm -e gdb -ex run cmake-build-debug/cppprojct 
