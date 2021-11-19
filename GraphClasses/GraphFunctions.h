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


template<typename DIFFEQ, typename SOLVER>
void contribute_to_integration(ReferenceContainer &REF, GeneralSolver<DIFFEQ,SOLVER> &solver, int MYPROC){
    PRINTF_DBG("Starting the contribution to int\n");std::cout<< std::flush;
    // Define timeout utilities
    const int DT= 5;
    int totlaps=0;
    // Build a graph vertex iterator
    auto vs = vertices(*REF.p_g);
    auto start = vs.first;
    auto end = vs.second;
    auto v = start;

    // Define relevant variables
    bool keepGoing = true;
    bool is_in = false;
    long * i;
    int remaining;
    long ix, currentuix;
    unsigned long uix;
    bool wait = false;

#pragma omp critical
    {
        remaining = *(REF.p_PENDING_INT);
    }

    // Main loop
    while (keepGoing) {
        // Fetch the data
#pragma omp critical
        {
            if (!REF.p_READY_FOR_INTEGRATION->second.empty()){
                ix = REF.p_READY_FOR_INTEGRATION->first.front();
                uix = REF.p_READY_FOR_INTEGRATION->second.front();
                REF.p_READY_FOR_INTEGRATION->first.pop();
                REF.p_READY_FOR_INTEGRATION->second.pop();
            } else {
                wait = true;
            }
        }

        if (wait) {
            mssleep(DT);
            PRINTF_DBG("(NOTHING TO INT)\n");std::cout<< std::flush;
        } else {
            ++totlaps;
            // Locate the graph vertex :-)
            PRINTF_DBG("starting to locate this ix: %lu\n",uix);
            v = start;
            while (!is_in) {
                currentuix = (unsigned long) get(get(boost::vertex_index, *(REF.p_g)), *v);
                if (currentuix == uix){
                    is_in = true;
                } else {
                    v++;
                }
                if (v == end) exit(1);
            }

            // Join the data as required by the integrator
            std::vector<double> neighborValues;
            std::vector<double> edgeValues;
            neighborValues.resize((*REF.p_IntHelper)[ix].ResultsPendProcess.size());
            edgeValues.resize((*REF.p_IntHelper)[ix].ResultsPendProcess.size());
            int _it = 0;
            for (auto const& it : (*REF.p_IntHelper)[ix].ResultsPendProcess) {
                neighborValues[_it] = std::get<0>(it);
                edgeValues[_it] = std::get<1>(it);
                _it++;
            }

            // Write the new value in the temporal register
            PRINTF_DBG("Integrating: temporal register was: %f", (*REF.p_g)[*v].temporal_register);
            double result = solver.evolve((*REF.p_IntHelper)[ix].centralValue,
                                          (*REF.p_IntHelper)[ix].centralParams,
                                          neighborValues,
                                          edgeValues);
            PRINTF_DBG("Integrating: temporal result was: %f", result);
            (*REF.p_g)[*v].temporal_register = result;
            PRINTF_DBG("Integrating: now it is: %f", (*REF.p_g)[*v].temporal_register);

            // Decrease the global value of pending integration
//#pragma omp critical
//{
//            *(REF.p_PENDING_INT)--;
//}
        }
        wait = false;

        // Communicate with other threads to see if we need to keep going
#pragma omp critical
{
        remaining = REF.p_READY_FOR_INTEGRATION->first.size();
}

        PRINTF_DBG("There are %d remaining vals to integrate...\n", remaining);std::cout<<std::flush;

        // Recompute the status of 'KeepGoing'
        keepGoing = (remaining != 0);
    }
    PRINTF_DBG("In total this thread of PROC %d performed %d laps!\n",MYPROC, totlaps);std::cout<<std::flush;
}


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


    std::pair<std::queue<long>, std::queue<unsigned long>> CHECKED;
    std::pair<std::queue<long>, std::queue<unsigned long>> READY_FOR_INTEGRATION;

    NT = ComHelper.NUM_THREADS;
    auto vs = vertices(g);
    double temporalResult;
    bool is_unclaimed = true;
    int TOT = 1;
    int active_responders=0, finalized_responders=0;
    bool keep_responding = true;
    ReferenceContainer REF(ParHelper,
                           ComHelper,
                           g,
                           CHECKED,
                           READY_FOR_INTEGRATION,
                           IntHelper,
                           TOT,
                           PENDING_INT,
                           MapHelper,
                           keep_responding);

    const int MAX_SUBTHR = 1;
    const int TIMETOL = 20;
    //ComHelper.WORLD_SIZE[0] / BATCH + 1;
    const int DT = 1;


    int request_performers=0;
    int request_performers_ended=0;

    // for tests
    std::set<int> Working, Finished;



#pragma omp parallel firstprivate(NVtot, vs, NT, N_total_nodes, REF, MAX_SUBTHR, TIMETOL, DT) // Node iterators have random access ;-)
    {
        bool am_i_first;
        int SplitCoef = omp_get_num_threads()/2; // 3/2
        if (SplitCoef < 2){SplitCoef = 2;}
        OpenMPHelper OmpHelper(NVtot, SplitCoef);


        if (OmpHelper.MY_THREAD_n < SplitCoef){

            if (OmpHelper.MY_THREAD_n % (MAX_SUBTHR + 1) == 0) {

                bool atomic_bool;
#pragma omp atomic read
                atomic_bool = keep_responding;

#pragma omp atomic update
                ++active_responders;
                    while (atomic_bool) { // as long as we keep processing our own,
                    //                    we mantain at least one dispatcher alive :-)
                        answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                        answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
#pragma omp atomic read
                        atomic_bool = keep_responding;
                        mssleep(DT);
                    }
#pragma omp atomic update
                ++finalized_responders;

                if (!atomic_bool) PRINTF_DBG("OVER!  :O\n");
                PRINTF_DBG(" I am thread %d (MESSAGE ANSWERER) and I have finished doing my job ;-)\n",omp_get_thread_num());

            } else {
#pragma omp atomic update
                request_performers++;
                perform_requests<DT, TIMETOL, BATCH>(NVtot, REF, N_total_nodes, OmpHelper);
#pragma omp atomic update
                request_performers_ended++;
            }
        } else {
            unsigned long NLocals, NInedges, M, rank, NOwned;
            long i=-1;
            unsigned long ui =0;

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
                if (get(MapHelper.NodeOwner,*v) != REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n]){
                    // we are treating a node which we dont own, please ignore it
                    PRINTF_DBG("\n\n\nIGNORING NODE THAT IS NOT OURS!\n\n\n\n");std::cout<<std::flush;
                    //++v;
                    exit(1);
                } else {
                    ui = ((unsigned long) get(get(boost::vertex_index, g), *v));
                }
                auto neighbors = boost::adjacent_vertices(*v, g);
                auto in_edges = boost::in_edges(*v, g);
#pragma omp parallel firstprivate(i, v, M, neighbors, NLocals, NInedges, NOwned, N_total_nodes, solver)
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
                        READY_FOR_INTEGRATION.first.push(i);
                        READY_FOR_INTEGRATION.second.push(ui);
                    }
#pragma omp atomic update
                    ++TOT;
                    PRINTF_DBG("Increased by one the number of total vertex\n");std::cout<<std::flush;
                } else {
#pragma omp critical
                    {
                        CHECKED.first.push(i); // Adding the index to the list of checked indexes ;-)
                    }
#pragma omp critical
                    {
                        CHECKED.second.push(ui); // Adding the index to the list of checked indexes ;-)
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

            // FIRST are we over: checking that all nodes have been processed
            while (!are_we_over){
                // Re-check if we are over
#pragma omp atomic read
                atomical_int = TOT;
                are_we_over = (atomical_int >= NVtot); // TODO: same here ;-).
                ++notreadyyet;
                mssleep(DT);
                PRINTF_DBG("not ready yet!\n");
            }

            // SECOND are we over: checking that all perform requesters have exited
            int aux1,aux2;
#pragma omp atomic read
            aux1 = request_performers;
#pragma omp atomic read
            aux2 = request_performers_ended;
            are_we_over = (aux1 == aux2);
            while (!are_we_over){
#pragma omp atomic read
                aux2 = request_performers_ended;
                are_we_over = (aux1 == aux2);
                mssleep(DT);
    //            printf("Padded lap in order to graciously terminate local performance_requests ;-)\n");
            }

            PRINTF_DBG("\n\n\n\n\n\n\WE WERE OVERRR NVtot and TOT are: %d & %d\n\n\n\n", NVtot, TOT);
            PRINTF_DBG("Before being ready, we waited for %d laps!\n", notreadyyet);
            std::cout << std::flush;


            int worldsize = ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n];
            int worldrank = ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n];

            PRINTF_DBG("\nabout to 1st barrier says P:%d\n", worldrank);std::cout<<std::flush;
            MPI_Barrier(MPI_COMM_WORLD);

            PRINTF_DBG("The (FIRST) cherry of the cake ;-)\n");std::cout<<std::flush;
#pragma omp atomic write
            keep_responding = false;

            bool readme;
#pragma omp atomic read
            readme = keep_responding;
            PRINTF_DBG("keep_responding was now set as %d\n",readme);std::cout<<std::flush;

            int auxy1, auxy2;
#pragma omp atomic read
            auxy1 = active_responders;
#pragma omp atomic read
            auxy2 = finalized_responders;
            int jh=0,jmax=1000;
            while (auxy1 > auxy2){
                mssleep(DT);
#pragma omp atomic read
                auxy2 = finalized_responders;
#pragma omp atomic read
                auxy1 = active_responders;
                PRINTF_DBG("active (total) are %d and finalized are %d\n", auxy1, auxy2);
                jh++;
                if (jh>jmax) PRINTF_DBG("active (total) are %d and finalized are %d\n", auxy1, auxy2);
            }
            PRINTF_DBG("\nabout to 2nd barrier says P:%d\n", worldrank);std::cout<<std::flush;
            MPI_Barrier(MPI_COMM_WORLD);
            PRINTF_DBG("\nEnded barrier says P:%d\n", worldrank);std::cout<<std::flush;

            PRINTF_DBG("The (SECOND) cherry of the cake (-;\n");std::cout<<std::flush;

        } else if (OmpHelper.MY_THREAD_n % 2 == 0) {
            int atomic_bool;
            bool later_mark_finalized = false;
            if (atomic_bool){
#pragma omp atomic update
                active_responders ++;
                later_mark_finalized = true;
                }
            while (atomic_bool){
            //for (int k=0;k<50;++k){
                //printf("A for worker that is finally responding messages says that we are not over with TOT\n");
                // 1) Answer some messages
                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);

                // 2) Re-check if we are over
#pragma omp atomic read
                atomic_bool = keep_responding;
                PRINTF_DBG("so far I keep responding U.u cuz keep responding was %d\n", atomic_bool);
            }
            if (later_mark_finalized) {
#pragma omp atomic update
                finalized_responders++;
            }
        }

        // Integration section
        PRINTF_DBG("starting to contribute\n");std::cout<<std::flush;
        int auxv1,auxv2;
#pragma omp atomic read
        auxv1 = request_performers;
#pragma omp atomic read
        auxv2 = request_performers_ended;
        bool keep_integrating = (auxv1 > auxv2);
        while (keep_integrating) {
            contribute_to_integration(REF, solver,
                                      REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n]);
#pragma omp atomic read
            auxv1 = request_performers;
#pragma omp atomic read
            auxv2 = request_performers_ended;
            keep_integrating = (auxv1 > auxv2);
            mssleep(DT);
        }
        contribute_to_integration(REF, solver,
                                  REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n]);
        PRINTF_DBG("ending to contribute\n");std::cout<<std::flush;


} // end of the parallel construct

    PRINTF_DBG("starting to swap register\n");std::cout<<std::flush;
    // Swap temporal and main registers
    register_to_value(g);

    PRINTF_DBG("exited");
    PRINTF_DBG("About to synchronize");std::cout<<std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    PRINTF_DBG("Done");
    printf("-^-^-H-E-A-R-T-B-E-A-T-^-^-");
    PRINTF_DBG("\n\n\n\\n\n\n\n\n\n");std::cout<<std::flush;
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H