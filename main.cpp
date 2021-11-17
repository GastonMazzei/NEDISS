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
    check_nested_status();
    int tagMaxFlag;
    int *tag_ub;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &tagMaxFlag);

    // Testing Section ;-)
    const int NNodes = std::stoi(std::getenv("NNODES")); //12345;200;
    certify_tagMax_compliant(tagMaxFlag, *tag_ub, NNodes, size, VERTEXVAL_REQUEST_FLAG);
    const int NRUNS = 500;
    const int BATCH = 1; // Keep batch at 1: increasing it enables sending messages of size bigger than 1, but it is not fully implemented in CommunicationFunctions.h.
    const double ErdosRenyiProba = 0.5;
    const int TESTN = std::stoi(std::getenv("TEST"));

    if (TESTN == 0){
        // TEST INITIALIZATION OF ALL NETWORKS INCLUDING THEIR ATTRIBUTES
        graph_tests_init(SEED, NNodes, ErdosRenyiProba);
    } else if (TESTN == 1) {
        // TEST SINGLESTEP EVOLUTION 'A COUPLE' OF TIMES
        graph_tests_singlestep_evolution<BATCH>(SEED,
                                            NNodes,
                                            ErdosRenyiProba);
    } else if (TESTN == 2){
        // TEST SINGLESTEP EVOLUTION SEVERAL THOUSAND TIMES
        central_test_long_singlestep_run<BATCH>(SEED,
                                            NNodes,
                                            ErdosRenyiProba,
                                            NRUNS);

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