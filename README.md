
<center><h1><b>V 1.0.0 is ready!</b></h1> <i>(Decent documentation coming soon!)</i></center>

---


<center><i>How to run</i></center>
<i>Just configure the simulation via `simulation.conf.sh` and execute the desired script</i>

1) `run.sh` runs the program in whichever mode was specified in the configuration file
2) `produce_graphical_simulation.sh` runs the program in simulation mode and produces a graphical simulation
4) `run_debugxterm.sh` runs with gdb and displays each processor in a different console
5) `run_pipe.sh` runs and pipes the output to a txt file
6) `run_valgrind.sh` runs and produces a valgrind report
7) `run_time.sh` runs nonverbose and computes the execution time

The configuration file allows the following to be specified:

* EQNUMBER: `{0,1}` indexes one of the ODEs described below. Currently `0` for 1D linear and `1` for Noiseless Kuramoto.

* SOLVER: `{0,1}` indexes the solver. `0` for the Euler method and `1` for Runge Kutta.

* SOLVERDEPTH: `{1,2,3,4}` is the depth used in Euler's method. A depth of i with i>1 requires the (i-1)-th derivative of the Equation's field to be implemented.

* K1, K2, K3, K4: `[0,1]` are the weights to apply to each of the Runge Kutta terms.

* TOPOLOGY: `{0,1,2,3}` indexes the topology. Currently `0` is a ring, `1` is a Clique, `2` is Erdos-Renyi, and `3` is Small-World. The latter currently has its probability p fixed at compilation. 

* NRUNS: `0-inf` defines the number of iterations to perform for the test that can be accessed through the flag TEST=2.

* NNODES: `P*NT*2 - inf` defines the number of total nodes to use. A pathological lower bound can be P*NT*2, which means two times the number of processors times the number of threads per processor.

* TEST: `{-1,1,2,3}` is the running mode, see below.

* SEED: `any int` is the seed to use during the entire simulation.

* kneigh: `0-inf` defines the number of neighbors to which connect if building a Small-World graph.

* A,B: `[0,1]` define the proba to initialize the ScaleFree graph.

* J,WMIN,WMAX `double` define the (constant) coupling strength for kuramoto, and the min & max bounds for uniformly sampling the natural frequency of each node.

* proba: `[0,1]` defines the probability used in the relevant constructors, currently Erdos-Renyi and Small-World.

* OMP_THREAD_LIMIT: `3-inf` is the maximum number of threads allowed per processor. At least three are required by our current division of tasks which is thread-id dependant.

* OMP_NESTED: `true` is a flag for OpenMP and it's required to be `true` .

* np is the number of processors defined  to mpirun, where the flag `-oversubscribe` allows to use more than the actual number of threads when ran locally in a laptop.

* SAMPLING_FREQ: `{1,...,NRUNS}` is the number of iterations to wait between recording the state of the system in a .dot file that can be later used for both graphical plotting and simple value-access. 

<center><i>Equations & Solvers</i></center>

| Equation | Euler order 1 | Euler up to O(4) | Generalized Runge Kutta|
| --- | --- | --- | --- |
| 1D linear | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| Noiseless <a href='https://en.wikipedia.org/wiki/Kuramoto_model'>Kuramoto</a> | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| N-D linear | :x: | :x: | :x: |

<center><i>Topology & Initialization</i></center>

| Graph | All nodes 1 value | Nodes different values | All edges 1 value | Edges different values | 
| --- | --- | --- | --- | --- |
| Ring(N) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :hourglass: |
| Clique(N) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :hourglass: |
| Erdos Renyi(N,p) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :hourglass: |
| ScaleFree(N,A,B) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :hourglass: |
| Small World(N,k,p) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :x: |
| Grid(N,M) | :hourglass: | :hourglass: | :hourglass: | :hourglass: |
| Torus(N,M) | :x: | :x: | :x: | :x: |
| hypercube(N) | :x: | :x: | :x: | :x: |


<center><i>Working modes</i></center>

-1 - Produces a graphical simulation

0 - Tests initialization

1 - Tests single-step evolution

2 - Tests results after several runs


<center><i>How to contribute</i></center>
<i>Open a pull request if you have any idea; at the present time (beggining Dec. 2021) the main focus of attention should be adding random walks to measure diffusion. Also a way to compute Synchronization time in-house (i.e. without producing graphical information) could be useful. Please contact </i>


<center><i>Other Useful Tools</i></center>
 
- compiling with the `VERBOSE` macro defined produces a loot of output in each run.
- `CONSIDERATIONS.txt` contains special considerations, e.g. how the MPI implementation can affect us.
- `Guideline-DiffEqsSolvers[ES].txt` is a proto explanation of how to build a Differential Equation for users. It says 'ES' because it is in Spanish! Ouch.
