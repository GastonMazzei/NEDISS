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
#include "Simulation/SimulationPipe.h"


using namespace std;

// FLAGGED FOR KILLING
//using namespace Eigen;
// Eigen is supposedly optimized for "using the entire processor",
// thus hyperthreading is counterproductive
// source: https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
//#define HYPERTHREADING 0      // 1 if hyperthreading is on, 0 otherwise




int main(int argc, char** argv)
{

    // Start the MPI Session, capture common variables
    // and enforce MPI Implementation constraints
    int rank, size;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        error_report("[error] The MPI did not provide the requested threading behaviour.");
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const unsigned long NNodes = std::stoi(std::getenv("NNODES"));
    static int OMP_THREAD_LIMIT = std::stoi(std::getenv("OMP_THREAD_LIMIT"));
    unsigned int SEED = std::stoi(std::getenv("SEED"));
    const int TESTN = std::stoi(std::getenv("TEST"));
    const int TOPOLOGY = std::stoi(std::getenv("TOPOLOGY"));
    SolverConfig solver_config;
    int tagMaxFlag;
    int *tag_ub;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &tagMaxFlag);
    certify_tagMax_compliant(tagMaxFlag, *tag_ub, NNodes, size, VERTEXVAL_REQUEST_FLAG);


    // Enforce OpenMP constraints
    if (OMP_THREAD_LIMIT<3){error_report(min_threads);};
    check_nested_status();

    // Processor 0 receives a welcome message in it's console :-)
    print_init(rank);
    print_warnings(rank);

    // Define the only hyperparameter
    // todo: communications.h appears to not support BATCH != 1
    const int BATCH = 1;

    // Core section :-)
    if (TESTN == 0){
        graph_tests_init(TOPOLOGY, SEED, NNodes);
    } else if (TESTN == 1) {
        int EQNUMBER = std::stoi(std::getenv("EQNUMBER"));
        graph_tests_singlestep_evolution<BATCH>(SEED,
                                            NNodes,
                                            solver_config, TOPOLOGY, EQNUMBER);
    } else if (TESTN == 2){
        const int NRUNS = std::stoi(std::getenv("NRUNS"));
        int EQNUMBER = std::stoi(std::getenv("EQNUMBER"));
        central_test_long_singlestep_run<BATCH>(SEED,
                                                NNodes,
                                                NRUNS, solver_config, TOPOLOGY, EQNUMBER);
    } else if (TESTN == -1){
        const int NRUNS = std::stoi(std::getenv("NRUNS"));
        int EQNUMBER = std::stoi(std::getenv("EQNUMBER"));
        simulator<BATCH>(SEED,
                            NNodes,
                            NRUNS, solver_config, TOPOLOGY, EQNUMBER);
    }


    // Report the end of the script
    if (rank == 0) {
        cout << "End of script has been reached" << endl;
    }

    // End MPI Session
    int exit_status = MPI_Finalize();
    if (exit_status != 0){
        cout << "MPI_Finalize() yielded " << exit_status << endl;
    }

    // End gracefully. This prevents valgrind from tagging memory blocks as "definetely lost"
    exit(0);
}