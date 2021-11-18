//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_GRAPHFUNCTIONS_H
#define CPPPROJCT_GRAPHFUNCTIONS_H
#include "../macros/macros.h"
#include "../Communication/CommunicationFunctions.h"
#include "../Solvers/GeneralSolver.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/memory_management.h"
#include "../Utils/msleep.h"
#include "../Utils/error.h"
#include "GeneralGraph.h"
#include "../Utils/HelperClasses.h"
#include <mpi.h>
#include <set>
#include <boost/mpi/environment.hpp>
#include "../Communication/CommunicationFunctions.h"



void contribute_to_integration(ReferenceContainer &REF);


template<int DT, int TIMETOL, int BATCH>
void answer_messages(ReferenceContainer &REF, int MYTHR);


template<int DT, int TIMETOL, int BATCH>
void answer_messages_edges(ReferenceContainer &REF,int MYTHR);


void sendReqForTest(int MYPROC, int i);

template <int DT, int TIMETOL, int BATCH>
void perform_requests(int NNodes,
                      ReferenceContainer REF,
                      unsigned long N,
                      OpenMPHelper &O);

// processors will iterate their "N=num_vertex(g)" nodes,
// and in each iteration they will either
//          (A)  (volume term)
//               iterate through their "NLocal=rank" locally-available neighbors
//               uploading the NLocal results into a container that can be used
//               to compute the central node's value using the solver evolver :-)
//          (B)  (surface term)
//               iterate through their "Nlocal<rank" locally-available neighbors,
//               while collecting "Nnonlocal" (what boost::graph::distributed call)
//               'vertex descriptors', 'edge descriptors', and the vertex's owner
//               in order to get involved in the following step which is requiring
//               a copy of this DynamicNode from the node's owner. Also, a new
//              iteration is performed over inward edges in order to collect the same
//              from them.
//
// It is in this context that we require keeping track of the surface (B) and volume (A) term's
// information using "N" (nodes in the process) containers that store a parallel cell object.

void register_to_value(Graph &g);

void destroyRequestWithoutCounter(MPI_Request &R);
void freeRequestWithoutCounter(MPI_Request &R);

// Fully declared template
template<typename DIFFEQ, typename SOLVER, int BATCH> // e.g. DIFFEQ = NoiselessKuramoto, SOLVER = EulerSolver<NoiselessKuramoto>
void single_evolution(Graph &g,
                      GeneralSolver<DIFFEQ,SOLVER> &solver,
                      CommunicationHelper &ComHelper,
                      ParallelHelper &ParHelper,
                      IntegrationHelper &IntHelper,
                      MappingHelper &MapHelper,
                      unsigned long N_total_nodes){

    // Roughly we are:
    // (1) gathering info from neighbors in parallel using MPI calls to other procs.
    // (2) solving the respective differential equation
    // (3) calling the proc synchronization

    unsigned long NVtot = boost::num_vertices(g), NT;
    int PENDING_INT = NVtot;


    std::queue<long> CHECKED, READY_FOR_INTEGRATION;

    NT = ComHelper.NUM_THREADS;
    auto vs = vertices(g);
    double temporalResult;
    bool is_unclaimed = true;
    int TOT = 1;
    ReferenceContainer REF(ParHelper,
                           ComHelper,
                           g,
                           CHECKED,
                           READY_FOR_INTEGRATION,
                           IntHelper,
                           TOT,
                           PENDING_INT,
                           MapHelper);

    const int MAX_SUBTHR = 1;
    const int TIMETOL = 20;
    //ComHelper.WORLD_SIZE[0] / BATCH + 1;
    const int DT = 1;
    bool keep_responding = true;

    // For testing :-0
    std::set<int> Working, Finished;
//    CHECKED.push(0);
//    CHECKED.push(1);
//    TOT = NVtot - 2;


#pragma omp parallel firstprivate(NVtot, vs, NT, N_total_nodes, REF, MAX_SUBTHR, TIMETOL, DT) // Node iterators have random access ;-)
{
    bool am_i_first;
    int SplitCoef = omp_get_num_threads()/2; // 3/2
    if (SplitCoef < 2){SplitCoef = 2;}
    OpenMPHelper OmpHelper(NVtot, SplitCoef);


    if (OmpHelper.MY_THREAD_n < SplitCoef){


        // SECTION OF NON-FOR-WORKERS :-)

        // EXPERIMENTAL:
//        if (OmpHelper.MY_THREAD_n == 0){
//            perform_requests<DT, TIMETOL, BATCH>(NVtot, REF, N_total_nodes,OmpHelper);
//        } else  if (OmpHelper.MY_THREAD_n == 1){
//            bool atomic_bool;
//#pragma omp atomic read
//            atomic_bool = keep_responding;
//            //for (int k=0; k<4; ++k) {
//            while (atomic_bool) { // as long as we keep processing our own,
//                //                    we mantain at least one dispatcher alive :-)
//                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
//                //answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
//#pragma omp atomic read
//                atomic_bool = keep_responding;
//            }
//        }


            // CLASSICAL
        if (OmpHelper.MY_THREAD_n % (MAX_SUBTHR + 1) == 0) {
            bool atomic_bool;
#pragma omp atomic read
            atomic_bool = keep_responding;
            //for (int k=0; k<4; ++k) {
            while (atomic_bool) { // as long as we keep processing our own,
                //                    we mantain at least one dispatcher alive :-)
                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                //answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
#pragma omp atomic read
                atomic_bool = keep_responding;
            }
            if (!atomic_bool) PRINTF_DBG("OVER!  :O\n");
            PRINTF_DBG(" I am thread %d (MESSAGE ANSWERER) and I have finished doing my job ;-)\n",omp_get_thread_num());
        } else { //if (OmpHelper.MY_THREAD_n == 1) {
            // DEBUG. change for just 'else' !!!
            perform_requests<DT, TIMETOL, BATCH>(NVtot, REF, N_total_nodes,OmpHelper);
            //printf("PERFORM_REQUEST HAS  ENDED! :-)");std::cout<<std::flush;
            //sendReqForTest(REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n], 0);
//#pragma omp atomic write
//            TOT = NVtot;

            PRINTF_DBG(" I am thread %d (PERFORM REQUESTER) and I have finished doing my job ;-)\n",omp_get_thread_num());
        }





    }
    else {
        unsigned long NLocals, NInedges, M, rank, NOwned;
        long i=-1;

        bool ready4int = false;
        i += OmpHelper.MY_OFFSET_n;
        for (auto v = vs.first + OmpHelper.MY_OFFSET_n;
             v != vs.first + OmpHelper.MY_OFFSET_n + OmpHelper.MY_LENGTH_n; ++v) {
            ++i;

            // Capture the central node's values
            IntHelper[i].centralValue = g[*v].value;
            IntHelper[i].centralParams = g[*v].params;
            IntHelper[i].build(g, *v, MapHelper, NOwned, rank, NLocals, M);
            if (NOwned == rank){
                ready4int = true;
            } else {
                ready4int = false;
            };

            auto neighbors = boost::adjacent_vertices(*v, g);
            auto in_edges = boost::in_edges(*v, g);
#pragma omp parallel firstprivate(i, v, M, neighbors, NLocals, NInedges, NOwned, N_total_nodes)
{
                OpenMPHelper OmpHelperN(NLocals, 0);
                for (auto n = neighbors.first + OmpHelperN.MY_OFFSET_n;
                     n != neighbors.first + OmpHelperN.MY_OFFSET_n + OmpHelperN.MY_LENGTH_n; ++n) {
                    auto local_e = edge(*v, *n, g).first;
                    auto local_v = *n;
                    if (get(MapHelper.Local, *n) == 1) { // Case (1): locally available elements ;-)
                        if (edge(*v, *n, g).second == 1) {
#pragma omp critical         // Just store the results directly in the results.
                            // we can probably afford being critical as there is
                           // MPI communication overhead in other threads.
{
                            IntHelper[i].ResultsPendProcess.emplace_back(g[*n].value,
                                                                         g[edge(*v, *n,
                                                                                g).first].value, // index to UID :-)
                                                                         ((unsigned long) get(
                                                                                 get(boost::vertex_index, g), *n)) *
                                                                         N_total_nodes);
}
                        } else { error_report("Push back mechanism for local nodes has failed"); };
                    } else { // Case (2.A): We 'see' this neighbor because it is connected
                        // via an edge that we own. Get the edge and record who is the other node's owner
#pragma omp critical
{
                        ParHelper.data[i].MissingA[OmpHelperN.MY_THREAD_n].emplace_back(g[local_e].value,
                                                                                        get(MapHelper.NodeOwner,
                                                                                            *n),
                                                                                        get(get(boost::vertex_index,
                                                                                                g), *n));
}
                    }
                }
                // Finally please contemplate the 3rd case:
                // we have in-edges that are not available
                // locally so we have neighbors that cant
                // be indexed by 'neighbors'
                OpenMPHelper OmpHelperE(M, 0, OmpHelperN.N_THREADS_n, OmpHelperN.MY_THREAD_n);
                int j = 0;
                long central_ix;
#pragma omp critical
{
                central_ix = get(get(boost::vertex_index,g), *v);
}
                for (auto e = in_edges.first; e != in_edges.second; ++e) {
                    if ((j >= OmpHelperE.MY_OFFSET_n) && (j < OmpHelperE.MY_OFFSET_n + OmpHelperE.MY_LENGTH_n)) {
                        auto local_e = *e;
                        auto local_v = boost::source(*e, g);
#pragma omp critical
{
                        ParHelper.data[i].MissingB[OmpHelperE.MY_THREAD_n].emplace_back(
                                                                                        (double) central_ix,
                                                                                        get(MapHelper.NodeOwner,
                                                                                            local_v),
                                                                                        get(get(boost::vertex_index,
                                                                                                g), local_v));
}
                    }
                    ++j;
                }
} // end of the nested parallelism! :-)


            if (ready4int){
#pragma omp critical
{
                READY_FOR_INTEGRATION.push(i);
}
#pragma omp atomic update
                ++TOT;
                PRINTF_DBG("Increased by one the number of total vertex\n");std::cout<<std::flush;
            } else {
#pragma omp critical
{
                CHECKED.push(i); // Adding the index to the list of checked indexes ;-)
}
                PRINTF_DBG("Added ix %d to CHECKED\n", i);std::cout<<std::flush;
            }
        }
        PRINTF_DBG(" I am thread %d (FOR WORKER) and I have finished doing my job ;-)\n",omp_get_thread_num());
        PRINTF_DBG("This thread is a for worker that has ended :-)\n");
    }

    //mssleep(15000); // keep guys waiting so we focus on MPI debugging ;_)

//    // BYPASS for debugging**
//#pragma omp atomic write
//    TOT = NVtot;
    // **********************

    //printf("I am a thread that  has  arrives to the end of the script :-)\n", ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n]);
    // The previous code is a race to this point :-) if you didnt arrived first move on lol
        PRINTF_DBG("about to enter 'omp critical unclaimed'\n"); std::cout << std::flush;
#pragma omp critical
{
            am_i_first = is_unclaimed;
            if (is_unclaimed){
                is_unclaimed = false;
                PRINTF_DBG("Claimed by thread %ld of processor %d\n", OmpHelper.MY_THREAD_n, ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n]);
            }
}

    if (am_i_first) {
        // Prepare vars
        int atomical_int;
        bool are_we_over = false;

        // While our Process is not over, I am the official responder :-)
#pragma omp atomic read
        atomical_int = TOT;
        PRINTF_DBG("Already checking if we are over, proc %d\n", ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n]);
        are_we_over = (atomical_int >= NVtot); // TODO: bug. it should be exactly equal (looks at perform_requests)
        int notreadyyet=0;
        while (!are_we_over){
            // 1) Answer some messages
            answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
            //answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
            // 2) Re-check if we are over
#pragma omp atomic read
            atomical_int = TOT;
            are_we_over = (atomical_int >= NVtot); // TODO: same here ;-).
            ++notreadyyet;
        }
        PRINTF_DBG("\n\n\n\n\n\n\WE WERE OVERRR NVtot and TOT are: %d & %d\n\n\n\n", NVtot, TOT);
        PRINTF_DBG("Before being ready, we waited for %d laps!\n", notreadyyet);
        std::cout << std::flush;

        // Once everyone in our Process is over, I can participate in the asynchronous
        // All2All while I keep responding.
        int counter=0;
        const int MAX_COUNTER=5;
        bool is_everyone_over = false;
        bool is_everyone_over_doubleChecked = false;
        int all2all_status = 0;
        int CHECKVAL;
        MPI_Request my_request[MAX_COUNTER]; // 5 requests
        int we_are_over[MAX_COUNTER][ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n]];
        int recv_status[MAX_COUNTER][ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n]];
        int return_status_all2all = 1;
        int return_status_getstatus = 1;
        std::set<int> ready;

        // Print for debug
        PRINTF_DBG("Already checking if everyone is over, proc %d, ready size is: %lu, here comes the ALL2ALL\n",
               ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n], ready.size());

        CHECKVAL = 1;
        while ((!is_everyone_over) || (!is_everyone_over_doubleChecked)) {

            // Currently not enforcing it as we dont increase the counter ;-)
            if (counter >= MAX_COUNTER){
                error_report("Max All2All instances exceeded");
            }
            if  (!is_everyone_over) {
                CHECKVAL = 1;
            } else  {
                CHECKVAL = 4;
                PRINTF_DBG("Now checkval is 4\n");
            }

            for (int i=0; i<ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n]; ++i) we_are_over[counter][i] = CHECKVAL;
            for (int i=0; i<ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n]; ++i) recv_status[counter][i] = 0;
            all2all_status = 0;
            ready = std::set<int>();
            ready.insert(ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n]);

            PRINTF_DBG("New all2all with checkval %d\n", CHECKVAL); // ampersand for recv_status after spotting RookieHPC's bug ;-)
            return_status_all2all = MPI_Ialltoall(&we_are_over[counter], 1, MPI_INT, &recv_status[counter], 1, MPI_INT, MPI_COMM_WORLD, &my_request[counter]);
            // Guarantee that the Iall2all is correctly generated
            while (return_status_all2all != 0) {
                PRINTF_DBG("return_status_all2all failed... answering messages before reattempting'\n"); std::cout << std::flush;
                // spend some time answering messages :-)
                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                //answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                // Reset the request
                PRINTF_DBG("about to free the necessary stuff in return status all2all\n"); std::cout << std::flush;
                freeRequestWithoutCounter(my_request[counter]);
                my_request[counter] = MPI_Request();
                // Try again
                PRINTF_DBG("about to reattemptsending return_status_all2all\n"); std::cout << std::flush;
                return_status_all2all = MPI_Ialltoall(&we_are_over[counter], 1, MPI_INT, &recv_status[counter], 1, MPI_INT, MPI_COMM_WORLD, &my_request[counter]);
            }

            // Guarantee that the status is correctly captured
            return_status_getstatus = MPI_Request_get_status(my_request[counter], &all2all_status, MPI_STATUS_IGNORE);
            while (return_status_getstatus != 0){
                PRINTF_DBG("return status getstatus failed... answering messages before reattempting'\n"); std::cout << std::flush;
                // spend some time answering messages :-)
                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                //answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);

                PRINTF_DBG("about to reattempt status getstatus\n"); std::cout << std::flush;
                return_status_getstatus = MPI_Request_get_status(my_request[counter], &all2all_status, MPI_STATUS_IGNORE);
            }

            // check if they have all been recieved :-)
            while (all2all_status != 1){
                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                //answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);

                MPI_Request_get_status(my_request[counter], &all2all_status, MPI_STATUS_IGNORE);
            }

            // if it has been recieved, update
            for (int pi = 0; pi < ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n]; pi++) {
                if (ready.count(pi) != 1) {
                    if (recv_status[counter][pi]==CHECKVAL) {
                        ready.insert(pi);
                    }
                }
            }
            if (ready.size() == ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n]) {
                if (is_everyone_over) {
                    is_everyone_over_doubleChecked = true;
                    PRINTF_DBG("P%d recieved consensus from all regarding that they are all over. Flags list is %d %d %d\n",
                           ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n],
                           recv_status[counter][0], recv_status[counter][1], recv_status[counter][2]);
                } else {
                    PRINTF_DBG("P%d says that everyone is over! flags list is %d %d %d\n",
                           ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n],
                           recv_status[counter][0], recv_status[counter][1], recv_status[counter][2]);
                    is_everyone_over = true;
                }
            } else {
                is_everyone_over = false;
                is_everyone_over_doubleChecked = false;
                // DEBUG MSG
                PRINTF_DBG("P%d says that not everyone has finished! flags list is %d %d %d\n",
                       ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n],
                       recv_status[counter][0], recv_status[counter][1], recv_status[counter][2]);
            }
            // Finally destroy request
            //destroyRequestWithoutCounter(my_request[counter]);
            freeRequestWithoutCounter(my_request[counter]);
            // Instead of increasing the counter, we keep using the same address
            my_request[counter] = MPI_Request();
            //++counter;
        }
        PRINTF_DBG("The cherry of the cake ;-)\n");std::cout<<std::flush;
#pragma omp atomic write
        keep_responding = false;

    } else  //if (OmpHelper.MY_THREAD_n % 2 == 0)
    {
        PRINTF_DBG("about to enter 'second branch of the end'\n"); std::cout << std::flush;
            // help answering messages ;-)
            bool atomic_bool;
#pragma omp atomic read
            atomic_bool = keep_responding;
            while (atomic_bool){
            //for (int k=0;k<50;++k){
                //printf("A for worker that is finally responding messages says that we are not over with TOT\n");
                // 1) Answer some messages
                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                //answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);

                // 2) Re-check if we are over
#pragma omp atomic read
                atomic_bool = keep_responding;
            }
        }
        contribute_to_integration(REF); // dont help answering messages ;-)
        // 2)
        //register_to_value(g);
} // end of the parallel construct
PRINTF_DBG("exited");
    PRINTF_DBG("About to synchronize");std::cout<<std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    printf("-^-Anti-Deadlock-Heartbeat-^-");std::cout<<std::flush;
    PRINTF_DBG("\n\n\n\\n\n\n\n\n\n");
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H