//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_GRAPHFUNCTIONS_H
#define CPPPROJCT_GRAPHFUNCTIONS_H
#include "../macros/macros.h"
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

template <int TIMETOL, int DT, int MAX_SUBTHR, int BATCH>
void GetAllMsgs(int NNodes,
                ReferenceContainer REF,
                unsigned long N,
                OpenMPHelper &O);

void contribute_to_integration(ReferenceContainer &REF);




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

    //void GetAllMsgs(int NNodes, CommunicationHelper &H, Graph &g, ParallelHelper &P, IntegrationHelper &I, std::queue<long> &C);
    unsigned long NVtot = boost::num_vertices(g), NT;
    int PENDING_INT = NVtot;


    std::queue<long> CHECKED, READY_FOR_INTEGRATION;

    NT = ComHelper.NUM_THREADS;
    auto vs = vertices(g);
    double temporalResult;

    int TOT = 1;
    ReferenceContainer REF(ParHelper,
                           ComHelper,
                           g,
                           CHECKED,
                           READY_FOR_INTEGRATION,
                           IntHelper,
                           TOT,
                           PENDING_INT);

    const int MAX_SUBTHR = 1;



#pragma omp parallel firstprivate(NVtot, vs, NT, N_total_nodes, REF, MAX_SUBTHR) // Node iterators have random access ;-)
{
    int SplitCoef = omp_get_num_threads()/2; // 3/2
    if (SplitCoef < 2){SplitCoef = 2;}
    OpenMPHelper OmpHelper(NVtot, SplitCoef);



    if (OmpHelper.MY_THREAD_n < SplitCoef){
                // The following is a templated function:
        // it requires <Number of timesteps, Delta time, Num Subthreads, BATCH>
        // we allow a delay of 5 failed attemps waiting 1ms and using only ONE helper (i.e. subthread)
        // the so-called 'BATCH' is how many requests are handled simultaneously (by each thread)
//        mssleep(5000);
        if (OmpHelper.MY_THREAD_n % (MAX_SUBTHR + 1) == 0) {
            GetAllMsgs<5, 1, MAX_SUBTHR, BATCH>(NVtot, REF, N_total_nodes, OmpHelper);
            GetAllMsgs<5, 1, MAX_SUBTHR, BATCH>(NVtot, REF, N_total_nodes, OmpHelper);
            // STRONG INTEGRATION OCCURS HERE! just call:
            if (VERBOSE) {
                // use PRINTF_DBG()
                std::cout << " I am thread " << omp_get_thread_num() <<
                          " and I have finished doing my job ;-)" << std::endl;
            }
            //contribute_to_integration(REF);
        }
    }
    else {
        //*************************************DEBUG!!!**************************************
//        if (OmpHelper.MY_THREAD_n == SplitCoef){
//#pragma omp critical
//{
//                CHECKED.push(0);
//                ParHelper.data[0].MissingA[0].emplace_back(4.,
//                                                                                1,
//                                                                                0);
//                CHECKED.push(1);
//                ParHelper.data[1].MissingA[0].emplace_back(5.,
//                                                                                2,
//                                                                                0);
//            }
//}
//        mssleep(50000); // DEBUG: watch what occurs in the MPI section :-0
        //*********************************************************************************
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

            // **********LITTLE CHEATING GOING ON HERE TO TEST COMMUNICATION STUFF
            ready4int = true;
            //**********************************************************************

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

            // FROZEN FOR DEBUGGING!
            if (ready4int){
#pragma omp critical
{
                READY_FOR_INTEGRATION.push(i);
}
#pragma omp atomic update
                --PENDING_INT;
#pragma omp atomic update
                ++TOT;
            } else {
#pragma omp critical
{
                CHECKED.push(i); // Adding the index to the list of checked indexes ;-)
}
            }
        } // end of the for :-)
        if (VERBOSE) {
            // use PRINTF_DBG()
            std::cout << " I am 'for' worker " << omp_get_thread_num() << " and I have completed my Job :-)" << std::endl;
        }
        // STRONG INTEGRATION OCCURS HERE! just call:
        contribute_to_integration(REF);
    }
} // end of the parallel construct

    // By now, there is only one thread and is not
    // dispatching other process' requests, so only
    // for assertion purposes we will produce:
    // wait until the global max of PENDING_INT is  exactly 0
    // use:
    if (PENDING_INT == 0) {

        // Define vars
        bool status = true;
        bool recv_status[ComHelper.WORLD_SIZE[0] - 1];
        int finalstate = MPI_Alltoall(&status,
                                      1,
                                      MPI_C_BOOL,
                                      &recv_status[0],
                                      ComHelper.WORLD_SIZE[0] - 1,
                                      MPI_C_BOOL,
                                      MPI_COMM_WORLD);
    } else
    {
        std::cout << " FAILED! check please" <<std::endl;
    }
    //std::cout << "Ended one lap!" << std::endl;
    //mssleep(5000);

    // swap (local) node's values with
    // the value at register
    register_to_value(g);


    // ************SYNCHRONIZATION DIRECTIVES: Only for Dynamic Networks
    // we are not modifying the topology so there is no need for synchronization
    //synchronize(g.process_group());
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H