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
#include <mpi.h>

void register_to_value(Graph &g);










// Fully declared template
template<typename DIFFEQ, typename SOLVER> // e.g. DIFFEQ = NoiselessKuramoto, SOLVER = EulerSolver<NoiselessKuramoto>
void single_evolution(Graph &g, GeneralSolver<DIFFEQ,SOLVER> &solver){
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

    // This could be useful:
    // https://www.boost.org/doc/libs/1_52_0/libs/graph/doc/adjacency_iterator.html
    //graph_traits<adjacency_list>::out_edge_iterator
    //graph_traits<adjacency_list>::adjacency_iterator
    //
    //boost::edge(u,v,g) returns pair<edge_descriptor, bool> where bool is if it exists

    // define mutable variables :-)
    typedef Graph::vertex_descriptor VD;
    typedef Graph::edge_descriptor ED;
    typedef std::pair<VD, ED> InfoVec_elem;
    typedef std::tuple<VD, ED, int> PartialInfoVec_elem;
    EdgeOwnerMap EdgeOwner = get(boost::edge_owner, g);
    OwnerMap NodeOwner = get(boost::vertex_owner, g);
    LocalVertexMap local = get(boost::vertex_local, g);
    GlobalVertexMap global = get(boost::vertex_global, g);
    auto vs = vertices(g);
    double centralValue, temporalResult;
    VD local_v; // dynamic pointers for vertex and edge descriptors
    ED local_e;  // (respectively)
    std::vector<double> centralParams;
    std::vector<double> edgeValues;
    std::vector<double> neighborValues;
    unsigned int ctrl_ownr, nghb_ownr, rank, MY_NUM = process_id(g.process_group());
    int i=-1, world_rank, world_size, isLocal, owner, NLocals, M;
    const unsigned long Nmax = boost::num_vertices(g);
    int lim_ProcessLocally[Nmax], counter_ProcessLocally = 0;
    int lim_MissingA[Nmax], counter_MissingA = 0;
    int lim_MissingB[Nmax], counter_MissingB = 0;
    InfoVec_elem ProcessLocally[Nmax][Nmax];
    PartialInfoVec_elem MissingA[Nmax][Nmax];
    PartialInfoVec_elem MissingB[Nmax][Nmax];
    //
    int DEBUG_FLAG = 0;
    //
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    boost::mpi::environment env;
    boost::mpi::communicator world;

    // Nous avons:
    // (VECTORS) centralParams, edgeValues, neighborValues
    // (propertyMap)    NodeOwner and EdgeOwner
    // (pair of vertex pointers) vs
    // (double) centralValue, temporalResult
    // (int) world_rank, world_size, MY_NUM, ctrl_ownr, nghb_ownr, i

    // Start the loop! :-)
    for (auto v = vs.first; v != vs.second; ++v) {
        ++i;

        // Get the central node's values
        centralValue = g[*v].value;
        centralParams = g[*v].params;

        //---------------USE-ME-IF-YOU-NEED-TO-DEBUG--;-)-----------------
        //std::cout << "So far so good, flag N " << DEBUG_FLAG << std::endl;
        //++DEBUG_FLAG;
        //----------------------------------------------------------------

        // Get the neighbor and edge's values
        auto neighbors = boost::adjacent_vertices(*v, g);
        rank = in_degree(*v, g) + out_degree(*v, g);
        NLocals = neighbors.second - neighbors.first;
        M = rank - NLocals;

        // THERE ARE 2 POSSIBLE SCENARIOS:
        // (1) all the data is available locally, which is Nlocals == rank
        // (2) there are M missing neighbors, where M = rank - NLocals,

        // Addressing (2):
        // there are 2 types of missing values.
        // 2.A) Nodes that are absent, but the edge object is available.
        //      This requires using that we know who is the owner and asking
        //      for the respective node.
        //      u=source(*ep.first,g);
        //      v=target(*ep.first,g);
        //      edge(vertex_descriptor u, vertex_descriptor v,
        //              const adjacency_list& g);
        //
        // 2.B) Nodes that are absent, and we only know about the incident edge.
        //      We must get the edge's owner and ask him for both the node and
        //      the edge's information.
        //      THIS DEFINES THE NEED FOR AN ASK_FOR_BOTH CALL

        if (M == 0) { // THE LEAST PROBABLE CASE! :-)
            // Skip directly to a Protocol to populate locally
//#pragma omp parallel for // Try me out! :-) Iterators are random access ones ;-)
            for (auto n = neighbors.first; n != neighbors.second; ++n) {
                local_e = edge(*v, *n, g).first;
                local_v = *n;
                ProcessLocally[i][counter_ProcessLocally].first = std::as_const(local_v);
                ProcessLocally[i][counter_ProcessLocally].second = std::as_const(local_e);
                counter_ProcessLocally++;
            }
            lim_ProcessLocally[i] = counter_ProcessLocally;
            counter_ProcessLocally = 0;
            counter_MissingA = 0;
            lim_MissingA[i] = counter_MissingA;
            counter_MissingB = 0;
            lim_MissingB[i] = counter_MissingB;
        } else { // THE MOST PROBABLE CASE! :-)
            // Iterate over neighbors
//#pragma omp parallel for
            for (auto n = neighbors.first; n != neighbors.second; ++n) {
                if (get(local, *n) == 1) { // Case (1): locally available elements ;-)
                    if (edge(*v, *n, g).second == 1) {
                        local_e = edge(*v, *n, g).first;
                        local_v = *n;
                        ProcessLocally[i][counter_ProcessLocally].first = std::as_const(local_v);
                        ProcessLocally[i][counter_ProcessLocally].second = std::as_const(local_e);
                        counter_ProcessLocally++;
                    } else { error_report("Push back mechanism for local nodes has failed"); };
                } else { // Case (2.A): We 'see' this neighbor because it is connected
                    // via an edge that we own. Get the edge and record who is the other node's owner
                    local_e = edge(*v, *n, g).first;
                    local_v = *n;
                    std::get<0>(MissingA[i][counter_MissingA]) = std::as_const(local_v);
                    std::get<1>(MissingA[i][counter_MissingA]) = std::as_const(local_e);
                    std::get<2>(MissingA[i][counter_MissingA]) = get(NodeOwner, *n);
                    counter_MissingA++;
                }
            }
            // Finally please contemplate the 3rd case:
            // we have in-edges that are not available
            // locally so we have neighbors that cant
            // be indexed by 'neighbors'
            auto in_edges = boost::in_edges(*v, g);
            for (auto e = in_edges.first; e != in_edges.second; ++e) {
                local_e = *e;
                local_v = boost::source(*e, g);
                std::get<0>(MissingB[i][counter_MissingB]) = std::as_const(local_v);
                std::get<1>(MissingB[i][counter_MissingB]) = std::as_const(local_e);
                std::get<2>(MissingB[i][counter_MissingB]) = get(EdgeOwner, *e);
                counter_MissingB++;
            }
            // Last but not least: clean the variables :-)
            lim_ProcessLocally[i] = counter_ProcessLocally;
            lim_MissingA[i] = counter_MissingA;
            lim_MissingB[i] = counter_MissingB;
            counter_ProcessLocally = 0;
            counter_MissingA = 0;
            counter_MissingB = 0;
        }
        adsync_barrier<70>();
        std::cout << "End of the neighbor iteration! starting new lap :-)" << std::endl;
        adsync_barrier<30>();
    }
    // At the end of the loop, we need to
    // set up a practical limit for future iterations:
    // now we know we only need ProcessLocally[0:lim_ProcessLocally[0], ..., 0:lim_ProcessLocally[P],0,..]
    // where P is counter_ProcessLocally = i (or maybe it's i-1?).
    counter_ProcessLocally = i;
    counter_MissingA = i;
    counter_MissingB = i;

    adsync_barrier<70>();
    std::cout << "End of the entire iteration! now we are ready to echange messages :-)" << std::endl;
    adsync_barrier<30>();


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




        // ******************HOW THE INTEGRATION CONTINUES  ;-)*****************
        //
        // Perform the evolution and store the result in the central
        // node's temporal register
//        g[*v].temporal_register = solver.evolve(centralValue,
//                                             centralParams,
//                                             neighborValues,
//                                             edgeValues);

        // Empty the vectors using a function from utils :-)
        //clear_vectors(centralParams, neighborValues, edgeValues);

    //}








    // We could block operations until all nodes have seen their neighbors ;-)
    // :-) there is no Queue so we must be strict.
    // Parallel BGL uses the BSP model, which we will enforce by replacing a barrier
    // with a synchronization directive.
    synchronize(g.process_group());
    //adsync_barrier<0>();


    // For all nodes,
    // central_value :=  temporal_register
    register_to_value(g);

    // Block operations until all nodes have been updated ;-)
    adsync_barrier<0>();
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H
