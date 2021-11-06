#include <iostream>

#include "mpi.h"
#include "Utils/timers.h"
#include "Tests/graph-test-init.h"
#include "Tests/graph-test-singlestep-evolution.h"
#include "Utils/global_standard_messages.h"
#include "Utils/print_init.h"
#include "Utils/print_warnings.h"
#include "Utils/parallel_sanitizer.h"

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
    print_init(rank);
    print_warnings(rank);
    static int OMP_THREAD_LIMIT = std::stoi(std::getenv("OMP_THREAD_LIMIT"));
    if (OMP_THREAD_LIMIT<3){error_report(min_threads);};
    check_nested_status();

    // Testing Section ;-)
    int N = 8;
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