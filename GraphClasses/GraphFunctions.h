//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_GRAPHFUNCTIONS_H
#define CPPPROJCT_GRAPHFUNCTIONS_H

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
template<typename DIFFEQ, typename SOLVER> // e.g. DIFFEQ = NoiselessKuramoto, SOLVER = EulerSolver<NoiselessKuramoto>
void single_evolution(Graph &g, GeneralSolver<DIFFEQ,SOLVER> &solver,
                      CommunicationHelper &ComHelper,
                      ParallelHelper &ParHelper,
                      IntegrationHelper &IntHelper,
                      MappingHelper &MapHelper){
    // There is an evolution operator that should compute the next
    // value and store it in the object's temporal register

    // said evolution operator should take (vector1, vector2, ...)
    //              where:
    //                  vector1:    <the value of the edges>
    //                  vector2:    <the value of the node at the other end of the edges>
    //                  vector3:    central node's parameters
    //                  value:      central node's value
    //                  ??MATRIX??: the parameter of each of the nodes, which should not be necessary
    //                              as it could be encoded in the edge.

    // TEST1: check that edges and nodes are correctly indexed, i.e.
    // (1)-1-(x)-2-(2) gets a vector1 and vector2 that is <1,2> and <1,2>

    // Instantiate the communication, parallel, integration and mapping helpers :-)
    //
    // !!! THEY SHOULD BE CONSUMED AS REFERENCES, NOT INSTANTIATED IN EACH CALL !!!
    //
    void GetAllMsgs(int NNodes, CommunicationHelper &H, Graph &g, ParallelHelper &P, IntegrationHelper &I, std::queue<long> &C);
    unsigned long NVtot = boost::num_vertices(g), NT;




    std::queue<long> CHECKED;
    NT = ComHelper.NUM_THREADS;
    auto vs = vertices(g);
    double temporalResult;

#pragma omp parallel firstprivate(NVtot, vs, NT)// Try me out! :-) Iterators are random access ones ;-)
{
    bool I_AM_NOT_MASTER = true; // this is mistrusting implementation of
                                // OpenMP standard where Master implies thread_id = 0 :-(
    long i=-1;
#pragma omp master
{
    I_AM_NOT_MASTER = false;
    GetAllMsgs(NVtot, ComHelper, g, ParHelper, IntHelper, CHECKED);
};
    if (I_AM_NOT_MASTER) {

        long NLocals, NInedges, M, rank;
        OpenMPHelper OmpHelper(NVtot, 1);


        i += OmpHelper.MY_OFFSET_n;
        for (auto v = vs.first + OmpHelper.MY_OFFSET_n;
                v != vs.first + OmpHelper.MY_OFFSET_n + OmpHelper.MY_LENGTH_n; ++v) {
            ++i;

            // Capture the central node's values
            IntHelper[i].centralValue = g[*v].value;
            IntHelper[i].centralParams = g[*v].params;
            // Get the neighbor and edge's values
            auto neighbors = boost::adjacent_vertices(*v, g);
            auto in_edges = boost::in_edges(*v, g);
            rank = in_degree(*v, g) + out_degree(*v, g);
            NLocals = neighbors.second - neighbors.first;
            M = rank - NLocals;

#pragma omp parallel firstprivate(i, M, neighbors, NLocals, NInedges)
{
            OpenMPHelper OmpHelperN(NLocals, 0);
            for (auto n = neighbors.first + OmpHelperN.MY_OFFSET_n;
                        n != neighbors.first + OmpHelperN.MY_OFFSET_n + OmpHelperN.MY_LENGTH_n; ++n) {
                auto local_e = edge(*v, *n, g).first;
                auto local_v = *n;
                if (get(MapHelper.Local, *n) == 1) { // Case (1): locally available elements ;-)
                    if (edge(*v, *n, g).second == 1) {
                        //std::make_tuple(std::as_const(local_v), std::as_const(local_e)));
                        ParHelper.data[i].ProcessLocally[OmpHelperN.MY_THREAD_n].emplace_back(local_v,
                                                                                    local_e,
                                                                                    get(get(boost::vertex_index, g), *v),
                                                                                    get(get(boost::vertex_index, g), *n));
                    } else { error_report("Push back mechanism for local nodes has failed"); };
                } else { // Case (2.A): We 'see' this neighbor because it is connected
                    // via an edge that we own. Get the edge and record who is the other node's owner
                    ParHelper.data[i].MissingA[OmpHelperN.MY_THREAD_n].emplace_back(local_v,
                                                                            local_e,
                                                                            get(get(boost::vertex_index, g), *v),
                                                                            get(get(boost::vertex_index, g), *n),
                                                                            get(MapHelper.NodeOwner, *n));
                }
            }
            // Finally please contemplate the 3rd case:
            // we have in-edges that are not available
            // locally so we have neighbors that cant
            // be indexed by 'neighbors'
            OpenMPHelper OmpHelperE(M, 0, OmpHelperN.N_THREADS_n, OmpHelperN.MY_THREAD_n);
            int j = 0;
            for (auto e = in_edges.first; e != in_edges.second; ++e) {
                if ((j>=OmpHelperE.MY_OFFSET_n) && (j<OmpHelperE.MY_OFFSET_n + OmpHelperE.MY_LENGTH_n)) {
                    auto local_e = *e;
                    auto local_v = boost::source(*e, g);
                    ParHelper.data[i].MissingB[OmpHelperE.MY_THREAD_n].emplace_back(local_v,
                                                                         local_e,
                                                                         get(get(boost::vertex_index, g), *v),
                                                                         get(get(boost::vertex_index, g), local_v),
                                                                         get(MapHelper.NodeOwner, local_v));
                }
                ++j;
            }
} // end of the nested parallelism! :-)
            std::cout << "End of the neighbor iteration! starting new lap :-)" << std::endl;
#pragma atomic
            CHECKED.push(i); // Adding the index to the list of checked indexes ;-)
        } // end of the for :-)
    } // end of the NON_MASTER_ONLY section :-)
} // end of the parallel construct

        // ******************HOW THE INTEGRATION CONTINUES  ;-)*****************
        //
        // Perform the evolution and store the result in the central
        // node's temporal register
//        g[*v].temporal_register = solver.evolve(centralValue,
//                                             centralParams,
//                                             neighborValues,
//                                             edgeValues);
//        clear_vectors(centralParams, neighborValues, edgeValues);





// **************************SOME COMMUNICATION WILL BE REQUIRED ;-)
//            if (world.rank() == 0) {
//                world.send(1, 0, std::string("Hello"));
//                std::string msg;
//                world.recv(1, 1, msg);
//                std::cout << msg << "!" << std::endl;
//            } else {
//                std::string msg;
//                world.recv(0, 0, msg);
//                std::cout << msg << ", ";
//                std::cout.flush();
//                world.send(0, 1, std::string("world"));
//            }

//            if ((world_rank == get(NodeOwner, *n)) && (world_rank != get(NodeOwner, *v))) {
////                MPI_Send(&g[*n].value, 1, MPI_DOUBLE, get(NodeOwner, *v), 0, MPI_COMM_WORLD);
//                world.send(get(NodeOwner, *v), 0, double(g[*n].value));
//            } else if ((world_rank == get(NodeOwner, *v)) && (world_rank != get(NodeOwner, *n))) {
//                std::cout << "I am process " << MY_NUM << " and I am trying to access a node I dont own "<< std::endl;
////                MPI_Recv(&temporalResult, 1, MPI_DOUBLE, get(NodeOwner, *n), 0, MPI_COMM_WORLD,
////                         MPI_STATUS_IGNORE);
//                world.recv(get(NodeOwner, *n), 0, temporalResult);
//                std::cout << "I recieved " << temporalResult << std::endl;
//                std::cout.flush();
//            }
    //adsync_barrier<20>();
    //}










    // We could block operations until all nodes have seen their neighbors ;-)
    // :-) there is no Queue so we must be strict.
    // Parallel BGL uses the BSP model, which we could enforce by replacing a barrier
    // with a synchronization directive.
    synchronize(g.process_group());


    // For all nodes,
    // central_value :=  temporal_register
    register_to_value(g);

    // Block operations until all nodes have been updated ;-)
    adsync_barrier<0>();
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H