

# ********************************** SIMULATION *******************************

# 0 Kuramoto, 1 LinearUncoupled
EQNUMBER=0;

# 0 init, 1 singlestep, 2 several steps, -1 print simulation
TEST=2;

SEED=21234;
SAMPLING_FREQ=15;
WMIN=1.3;
WMAX=0.3;
J=1.3;
NRUNS=4000;

# ********************************** SOLVER *******************************

# 0 EulerSolver, 1 RungeKutta
SOLVER=1;

# 1,2,3,4 for EulerSolver, or 1,2,4 for RungeKutta
SOLVERDEPTH=2;

# Coefficients for Runge Kutta
if [ $SOLVERDEPTH == 4 ]; then
 K1=0.1666666666666666;
 K2=0.33333333333333333;
 K3=0.333333333333333;
 K4=0.16666666666666;
elif [ $SOLVERDEPTH == 2 ]; then
 K1=0.500000000000000;
 K2=0.500000000000000;
 K3=0.000000000000000;
 K4=0.000000000000000;
else
 K1=1.000000000000000;
 K2=0.000000000000000;
 K3=0.000000000000000;
 K4=0.000000000000000;
fi

# ********************************** TOPOLOGY *******************************

# 0 Ring, 1 Clique, 2 ErdosRenyi, 3 SmallWorld
TOPOLOGY=3;
NNODES=500;
proba=0.15;
kneigh=50;

# ********************************** PERFORMANCE *******************************

OMP_THREAD_LIMIT=3;
OMP_NESTED=true;
NP=3;



