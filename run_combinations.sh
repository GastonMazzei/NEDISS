# Run several combinations to quickly verify that nothing has been broken


#
# Linear eq with Euler depth 4, Runge kutta (not implemented), Euler depth 1. Ring, Clique, ErdosRenyi were tried.
#

mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=0 -x SOLVERDEPTH=4 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=0 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct

echo "EXPECTED STATUS: OK"
sleep 1s

mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=1 -x SOLVERDEPTH=3 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=1 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct

echo "EXPECTED STATUS: placeholder"
sleep 1s

mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=0 -x SOLVERDEPTH=1 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=2 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct

echo "EXPECTED STATUS: OK"
sleep 1s



#
# Kuramoto eq with Euler depth 4, Runge kutta (not implemented), Euler depth 1. Ring, Clique, ErdosRenyi were tried.
#

mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=0 -x SOLVERDEPTH=4 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=0 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct

echo "EXPECTED STATUS: Error not implemented the derivatives"
sleep 1s

mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=1 -x SOLVERDEPTH=3 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=1 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct

echo "EXPECTED STATUS: placeholder"
sleep 1s

mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=0 -x SOLVERDEPTH=1 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=2 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct

echo "EXPECTED STATUS: OK"
sleep 1s


