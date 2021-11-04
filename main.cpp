#include <iostream>

#include "mpi.h"
#include "Utils/timers.h"
#include "Tests/graph-test-init.h"
#include "Tests/graph-test-singlestep-evolution.h"
#include <omp.h>
#include "Utils/typed_combinations.h"

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

int main(int argc, char** argv)
{

    // START:
    //
    int rank;
    unsigned int SEED = 12345;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Testing Section ;-)
    int N = 12;
    double p = 0.5;
    //graph_tests_init(SEED, N, p);
    graph_tests_singlestep_evolution(SEED, N, p);

    // END:
    //
    if (rank == 0) {
        cout << "End of script has been reached" << endl;
    }
    int exit_status = MPI_Finalize();
    return exit_status;
}