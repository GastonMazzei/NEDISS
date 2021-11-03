//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_GRAPHFUNCTIONS_H
#define CPPPROJCT_GRAPHFUNCTIONS_H

#include "../Solvers/GeneralSolver.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/memory_management.h"
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
    auto vs = vertices(g);
    double centralValue, temporalResult;
    typedef boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
    std::vector<double> centralParams;
    std::vector<double> edgeValues;
    std::vector<double> neighborValues;
    OwnerMap owner = get(boost::vertex_owner, g);

    unsigned int MY_NUM = process_id(g.process_group());
    unsigned int ctrl_ownr, nghb_ownr;
    int i;


    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    //std::cout  << "DEBUG!"<< std::endl;


//    auto v0 = vs.first;
//    auto neighbors = boost::adjacent_vertices(*v0, g);
//    ctrl_ownr = get(owner, *v0);
//    if (MY_NUM == nghb_ownr){
//        std::cout << "There are " << vs.second  - vs.first << " neighbors" << std::endl;
//    }
//    double results[4];
//    i=0;
//    for (auto n = neighbors.first; n != neighbors.second; ++n) {
//        nghb_ownr =  get(owner, *n);
////        if ((MY_NUM == nghb_ownr) && (MY_NUM == ctrl_ownr)) {
////            results[i] = g[*n].value;
////            std::cout << "First condition :-)" << std::endl;
////        } else {
////            if (MY_NUM == nghb_ownr) {
////                MPI_Send(&g[*n].value, 1, MPI_DOUBLE, ctrl_ownr, 0, MPI_COMM_WORLD);
////            } else if (MY_NUM == ctrl_ownr) {
////                MPI_Recv(&results[i], 1, MPI_DOUBLE, nghb_ownr, 0, MPI_COMM_WORLD,
////                         MPI_STATUS_IGNORE);
////                printf("Process %d received number %d from process %d\n",
////                       ctrl_ownr, results[i], nghb_ownr);
////            }
////            std::cout << "Second condition :-)" << std::endl;
////        }
////        if (MY_NUM == nghb_ownr){
////            std::cout << "Loop n " << i << std::endl;
////            std::cout << "Node owner is P: "<< nghb_ownr << std::endl;
////        }
//        std::cout << "I am process" << MY_NUM << " And I have noted one node" << std::endl;
//        ++i;
//    };




    for (auto v = vs.first; v != vs.second; ++v) {
        i = 0;
        // Get the central node's values
        centralValue = g[*v].value;
        centralParams = g[*v].params;
        ctrl_ownr = get(owner, *v);
        std::cout << "-MacroLap-" << std::endl;

        // Get the neighbor and edge's values
        auto neighbors = boost::adjacent_vertices(*v, g);
        for (auto n = neighbors.first; n != neighbors.second; ++n) {
            if (MY_NUM == get(owner, *n)) {
                neighborValues.push_back(g[*n].value);
                edgeValues.push_back(g[boost::edge(*v, *n, g).first].value);
            }
            std::cout << "I am process" << MY_NUM << " And I have noted one node. Am I its owner? " <<  get(owner, *n) << std::endl;
        }
        // Perform the evolution and store the result in the central
        // node's temporal register
//        g[*v].temporal_register = solver.evolve(centralValue,
//                                             centralParams,
//                                             neighborValues,
//                                             edgeValues);

        // Empty the vectors using a function from utils :-)
        //clear_vectors(centralParams, neighborValues, edgeValues);

    }

    // Block operations until all nodes have seen their neighbors ;-)
    adsync_barrier<0>();

    // For all nodes,
    // central_value :=  temporal_register
    register_to_value(g);

    // Block operations until all nodes have been updated ;-)
    adsync_barrier<0>();
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H
