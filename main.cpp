#include <iostream>

#include "mpi.h"
#include "Utils/timers.h"
#include "Utils/reproductibility.h"
#include "Tests/graph-test-init.h"
#include "Tests/graph-test-singlestep-evolution.h"
#include <omp.h>

//#include <string>

//#include  <cmath>
//#include <random>
//#include <algorithm>
//#include <iterator>
//#include <chrono>
//#include "Utils/error.h"
//#include <Eigen/Dense>
//#include <iostream>
//#include <math.h>
//#include <random>
//#include <vector>
//#include <iostream>

using namespace std;


//using namespace Eigen;
// Eigen is supposedly optimized for "using the entire processor",
// thus hyperthreading is counterproductive
// source: https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
//#define HYPERTHREADING 0      // 1 if hyperthreading is on, 0 otherwise

void graph_tests_init(unsigned int SEED){
    // ---CONSTRUCTOR AND INITIALIZATION DEBUGGING---

    // Clique Network
    reproductibility_lock(SEED);
    test_clique_graph_init(4);

    // Ring Network
    reproductibility_lock(SEED);
    test_ring_graph_init(4);

    // Erdos Renyi Network
    reproductibility_lock(SEED);
    test_erdosRenyi_graph_init(5, 0.5);
}


void graph_tests_singlestep_evolution(unsigned int SEED){
    // ---SINGLE TIMESTEP EVOLUTION DEBUGGING---

    // Ring Network
    reproductibility_lock(SEED);
    test_ring_graph_singlestep_evolution(4);
}

int main(int argc, char** argv)
{

    // START:
    //
    int rank;
    unsigned int SEED = 12345;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Testing Section ;-)
    //graph_tests_init(SEED);
    graph_tests_singlestep_evolution(SEED);

    // END:
    //
    if (rank == 0) {
        cout << "End of script has been reached" << endl;
    }
    int exit_status = MPI_Finalize();
    return exit_status;
}