#include <iostream>

#include "mpi.h"
#include "utils/timers.h"
#include "utils/reproductibility.h"
#include "tests/graph-test.h"
#include <omp.h>

//#include <string>

//#include  <cmath>
//#include <random>
//#include <algorithm>
//#include <iterator>
//#include <chrono>
//#include "utils/error.h"
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

void graph_tests(unsigned int SEED){

    // Clique Network
    reproductibility_lock(SEED);
    test_clique_graph(4);

    // Erdos Renyi Network
    //reproductibility_lock(SEED)
    //test_erdosRenyi_graph(5, 0.5);
}


int main(int argc, char** argv)
{

    // START:
    //
    int rank;
    unsigned int SEED = 12345;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    graph_tests(SEED);


    // END:
    //
    if (rank == 0) {
        cout << "End of script has been reached" << endl;
    }
    int exit_status = MPI_Finalize();
    return exit_status;
}