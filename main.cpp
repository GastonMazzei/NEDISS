#include "macros/macros.h"

#include <iostream>

#include "mpi.h"
#include "Utils/timers.h"
#include "Tests/graph-test-init.h"
#include "Tests/graph-test-singlestep-evolution.h"
#include "Utils/global_standard_messages.h"
#include "Utils/print_init.h"
#include "Utils/print_warnings.h"
#include "Utils/parallel_sanitizer.h"
#include "Tests/long-singlestep-run.h"
#include "Utils/certify_compliant.h"
#include "Simulation/SimulationExample.h"


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
    int rank, size;
    int provided;
    unsigned int SEED = std::stoi(std::getenv("SEED"));//12345;

    // MPI_Init(&argc, &argv); // Not this, we want a multithreaded MPI application :-)
    //
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        error_report("[error] The MPI did not provide the requested threading behaviour.");
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    print_init(rank);
    print_warnings(rank);
    static int OMP_THREAD_LIMIT = std::stoi(std::getenv("OMP_THREAD_LIMIT"));
    if (OMP_THREAD_LIMIT<3){error_report(min_threads);};
    //if (OMP_THREAD_LIMIT<3){error_report(min_threads);}; APPLY THIS ALSO TO 3 PROCESSORS MIN, because of debug prints that use [0][1][2] indexes ;-)
    check_nested_status();
    int tagMaxFlag;
    int *tag_ub;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &tagMaxFlag);

    // Testing Section ;-)
    const unsigned long NNodes = std::stoi(std::getenv("NNODES")); //12345;200;
    certify_tagMax_compliant(tagMaxFlag, *tag_ub, NNodes, size, VERTEXVAL_REQUEST_FLAG);
    const int BATCH = 1; // Keep batch at 1: increasing it enables sending messages of size bigger than 1, but it is not fully implemented in CommunicationFunctions.h.
    const int TESTN = std::stoi(std::getenv("TEST"));
    // Solvers: 0 is Euler, 1 is RungeKutta
    SolverConfig solver_config;
    // Topology: 0 is Ring, 1 is Clique, 2 is Erdos-Renyi, 3 is Small-World
    const int TOPOLOGY = std::stoi(std::getenv("TOPOLOGY"));

    if (TESTN == 0){
        // TEST INITIALIZATION OF ALL NETWORKS INCLUDING THEIR ATTRIBUTES
        graph_tests_init(TOPOLOGY, SEED, NNodes);
    } else if (TESTN == 1) {
        // TEST SINGLESTEP EVOLUTION 'A COUPLE' OF TIMES
        int EQNUMBER = std::stoi(std::getenv("EQNUMBER"));
        graph_tests_singlestep_evolution<BATCH>(SEED,
                                            NNodes,
                                            solver_config, TOPOLOGY, EQNUMBER);
    } else if (TESTN == 2){
        // TEST SINGLESTEP EVOLUTION SEVERAL THOUSAND TIMES
        const int NRUNS = std::stoi(std::getenv("NRUNS"));
        int EQNUMBER = std::stoi(std::getenv("EQNUMBER"));
        central_test_long_singlestep_run<BATCH>(SEED,
                                                NNodes,
                                                NRUNS, solver_config, TOPOLOGY, EQNUMBER);
    } else if (TESTN == -1){
        // TEST SINGLESTEP EVOLUTION SEVERAL THOUSAND TIMES
        const int NRUNS = std::stoi(std::getenv("NRUNS"));
        int EQNUMBER = std::stoi(std::getenv("EQNUMBER"));
        simulator<BATCH>(SEED,
                            NNodes,
                            NRUNS, solver_config, TOPOLOGY, EQNUMBER);
    }


    // END:
    //
    if (rank == 0) {
        cout << "End of script has been reached" << endl;
    }
    int exit_status = MPI_Finalize();

    exit(0); // prevents valgrind from tagging it as definetely lost :O
    //return exit_status;
}