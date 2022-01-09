
<center><h1><b>V 1.0.0</b></h1> 

 [![Build status](https://ci.appveyor.com/api/projects/status/8t72t8yb59kbdiuj?svg=true)](https://ci.appveyor.com/project/GastonMazzei/nediss)
[![Language grade: C++](https://img.shields.io/lgtm/grade/cpp/github/GastonMazzei/NEDISS)](https://lgtm.com/projects/g/GastonMazzei/NEDISS/context:cpp)

 
 
 <a href='https://gastonmazzei.github.io/NEDISS/annotated.html'><i>DOCS</i></a></center>
 
---


<center><b>How to run</b></center>

<i>Configure the simulation at `simulation.conf.sh` and execute the desired script. If it's a sweep, configuring `period_analysis_configuration.py` is also required. (TODO: video tutorial)</i>

1) `produce_graphical_simulation.sh` runs the program in simulation mode and produces a graphical simulation
2) `capture_period.sh` computes the exponential coefficients of the collapse into a locked state for the values specified in the configuration file.
3) `period_analysis.sh` wraps the capture of the period for a range of parameters indicated in the sweep configuration file, allows to mark outliers to exclude, and finally produces a chart showing the trend. 
4) TODO: compute the average return probability using Monte Carlo solver.


<center><b>Equations & Solvers</b></center>

| Equation | Euler order 1 | Euler up to O(4) | Generalized Runge Kutta | Solver Wrapper w/Noise models | Monte Carlo Evolution |
| --- | --- | --- | --- | --- | --- |
| 1D linear | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :x: | :bulb: | :x: |
| Noiseless <a href='https://en.wikipedia.org/wiki/Kuramoto_model'>Kuramoto</a> | :heavy_check_mark: | :x: | :heavy_check_mark: | :x: | :bulb: | :x: |
| N-D linear |:bulb: | :bulb: | :bulb: | :x: | :bulb: | :x: |
| <a href='http://www.scholarpedia.org/article/Pulse_coupled_oscillators'>Pulse-coupled</a> oscillators | :bulb: | :bulb: | :bulb: | :x:  | :bulb: | :x: |
| Fractional Brownian Motion | :x: | :x: | :x: | :x: | :bulb: |

|<b>Refs.</b>
| :heavy_check_mark: done and tested
| :x: not supported 
|:bulb: potential feature|


<center><b>Topologies & Initializations</b></center>

| Graph | All nodes 1 value | Nodes different random values | All edges 1 value | Edges different random values | Node and edges read from file | 
| --- | --- | --- | --- | --- | --- |
| Ring(N) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :bulb: | :bulb: |
| Clique(N) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :bulb: | :bulb: |
| Erdos Renyi(N,p) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :bulb: | :bulb: |
| ScaleFree(N,A,B) | :white_check_mark: | :white_check_mark: | :white_check_mark: | :bulb: | :bulb: |
| Small World(N,k,p) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :hourglass: | :hourglass: |
| Grid(N,M) | :bulb: | :bulb: | :bulb: | :bulb: | :bulb: |
| Torus(N,M) | :bulb: | :bulb: | :bulb: | :bulb: | :bulb: |
| hypercube(N) | :bulb: | :bulb: | :bulb: | :bulb: | :bulb: |

|<b>Refs.</b>
| :heavy_check_mark: done and tested
| :white_check_mark: done but untested 
|:hourglass: under construction 
|:bulb: potential feature|



<center><b>How to Configure</b></center>
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




<center><b>TODO: How to write plugins</b></center>

- Write your own initial condition

- Write your own equation

-  Write your own Graph Topology



<center><b>Other running modes</b></center>

- `run.sh` just runs the program in whichever mode was specified in the configuration file

- `run_debugxterm.sh` runs with gdb and displays each processor in a different console

- `run_pipe.sh` runs and pipes the output to a txt file

- `run_valgrind.sh` runs and produces a valgrind report

- `run_time.sh` runs nonverbose and computes the execution time


<center><b>Internal Working Modes (exposed to `run.sh`)</b></center>

-1 - Produces a simulation of length `NRUNS` and samples the node values at `SAMPLING_FREQ`

0 - Visually inspect node and edge values at initialization

1 - Visually inspect node and edge values after "a few" iterations (currently fixed at 3).

2 - Visually inspect node and edge values at both initialization and after several configurable (`NRUNS`) iterations. It also verifies the probability of errors, e.g. verification of 10k iterations with 100 threads and no segmentation fault.



<center><b>How to contribute</b></center>
<i>Please open a pull request if you have any idea; at the present time (beggining Jan. 2022) the main focus of attention should be adding random walks to measure diffusion. Do not hesitate to contact me. As a developer please consider compiling with the `VERBOSE` macro defined. The pattern for pushing minor changes without triggering an appveyor job is the commit title "[unsurveilled]".</i>
