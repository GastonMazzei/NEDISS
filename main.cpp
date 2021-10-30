#include <iostream>

#include "mpi.h"
#include "utils/timers.h"
#include "tests/graph-test.h"

//#include <string>
//#include <omp.h>
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


void graph_tests(){

    // Clique Network
    test_clique_graph(15);

    // Erdos Renyi Network
    test_erdosRenyi_graph(15, 0.3);
}


int main(int argc, char** argv)
{

    // START:
    //
    MPI_Init(&argc, &argv);

    graph_tests();

    // END:
    //
    cout << "End of script has been reached" << endl;
    int exit_status = MPI_Finalize();
    return exit_status;
}