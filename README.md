

A brief graphic explanation.

<img src="https://github.com/GastonMazzei/NEDISS/raw/main/front.png" width=800>


<center><i>How to run</i></center>

```mpirun -oversubscribe -x EQNUMBER=1 -x SOLVER=0 -x SOLVERDEPTH=3 -x K1=0.2 -x K2=0.4 -x K3=0.4 -x K4=0.2  -x TOPOLOGY=0 -x NRUNS=50 -x NNODES=25 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct```

* EQNUMBER: `{0,1}` indexes one of the ODEs described below. Currently `0` for 1D linear and `1` for Noiseless Kuramoto.

* SOLVER: `{0,1}` indexes the solver. `0` for the Euler method and `1` for Runge Kutta.

* SOLVERDEPTH: `{1,2,3,4}` is the depth used in Euler's method. A depth of i with i>1 requires the (i-1)-th derivative of the Equation's field to be implemented.

* K1, K2, K3, K4: `[0,1]` are the weights to apply to each of the Runge Kutta terms.

* TOPOLOGY: `{0,1,2}` indexes the topology. Currently `0` is a ring, `1` is a Clique, and `2` is Erdos-Renyi. The latter currently has its probability p fixed at compilation. 

* NRUNS: `0-inf` defines the number of iterations to perform for the test that can be accessed through the flag TEST=2.

* NNODES: `P*NT*2 - inf` defines the number of total nodes to use. A pathological lower bound can be P*NT*2, which means two times the number of processors times the number of threads per processor.

* TEST: `{0,1,2}` are the possible test modes to run as it is currently under development.

* SEED: `any int` is the seed to use during the entire simulation.

* OMP_THREAD_LIMIT: `3-inf` is the maximum number of threads allowed per processor. At least three are required by our current division of tasks which is thread-id dependant.

* OMP_NESTED: `true` is a flag for OpenMP and it's required to be `true` .

* np is the number of processors defined  to mpirun, where the flag `-oversubscribe` allows to use more than the actual number of threads when ran locally in a laptop.

* TODO: more inputs. (1) Erdos-Renyi proba, (2) Hyperparams to be optimized such as BATCH, DT.


<center><i>Equations & Solvers</i></center>

| Equation | Euler order 1 | Euler up to O(4) | Generalized Runge Kutta|
| --- | --- | --- | --- |
| 1D linear | :heavy_check_mark: | :heavy_check_mark: | :x: |
| Noiseless <a href='https://en.wikipedia.org/wiki/Kuramoto_model'>Kuramoto</a> | :heavy_check_mark: | :x: | :x: |
| N-D linear | :x: | :x: | :x: |

<center><i>Topology & Initialization</i></center>

| Graph | All nodes 1 value | Nodes different values | All edges 1 value | Edges different values | 
| --- | --- | --- | --- | --- |
| Ring(N) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :x: |
| Clique(N) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :x: |
| Erdos Renyi(N,p) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :black_square_button: |
| Grid(N,M) | :x: | :x: | :x: | :x: |
| Torus(N,M) | :x: | :x: | :x: | :x: |
| Small World(N,k,p) | :x: | :x: | :x: | :x: |
| hypercube(N) | :x: | :x: | :x: | :x: |


<center><i>Tests</i></center>

0 - Tests initialization
1 - Tests single-step evolution
2- Tests results after several runs

 
