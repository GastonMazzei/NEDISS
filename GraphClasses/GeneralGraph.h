//
// Created by m4zz31 on 29/10/21.
//

#ifndef CPPPROJCT_GENERALGRAPH_H
#define CPPPROJCT_GENERALGRAPH_H

#include  <iostream>

// Flagged to kill after the next test 01 jan 2022.
//#include <Eigen/Sparse>
#include <string>
#include "../Utils/error.h"

#include "../macros/macros.h"

#include <omp.h>
#include <vector>
#include <algorithm>



#include <boost/mpi.hpp>
#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/named_graph.hpp>
#include <boost/graph/parallel/process_group.hpp>
#include <boost/assert.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/parallel/caching_property_map.hpp>
#include <boost/graph/parallel/algorithm.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/parallel/process_group.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/parallel/local_property_map.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/named_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/range/adaptors.hpp>


struct DynamicNode {
    /*
    OBJECT DynamicNode are the elements we want to evolve over time
    it should have:
    (1) its value:
      - 1 dimensional for Kuramoto
      - TODO: N-dimensional (e.g. N=3 for Rossler Oscillators)
    (2) parameters:
      - an N-dimensional array with its own properties: characteristic frequency for example
    (3) temoral register:
      - 1 dimensional float to replace the value whenever the next update call arrives
    */

    double value = 0;
    double temporal_register = 0;
    std::vector<double> params;

    DynamicNode() = default;
    DynamicNode(double i): value(i) {};

    // Boost::graph requires the nodes to have a serialization attribute :-)
    template<typename Archiver>
    void serialize(Archiver& ar, const unsigned int) {
        ar & value & params ;
    }
};

struct DynamicEdge {
    /*
    Edges represent interactions between nodes
    The Kuramoto equation made each pair interact only once, thus a double
    value was deemed appropiate. The most general equation could require
    several interaction coefficients, e.g.
    df1/dx = coefficient1 * f1 * f2 + coefficient2 * f1 * f2^2
    todo: N-dimensional values for edges
    */
    double value = 1;

    DynamicEdge() = default;
    DynamicEdge(double i): value(i) {};

    // Boost::graph requires the nodes to have a serialization attribute :-)
    template<typename Archiver>
    void serialize(Archiver& ar, const unsigned int /*version*/) {
        ar & value;
    }
};

// flagged to kill after the next test 01 jan 2022.
//typedef DynamicNode DynamicNode;
//typedef DynamicEdge DynamicEdge;


typedef boost::adjacency_list<boost::vecS,
        boost::distributedS<boost::graph::distributed::mpi_process_group, boost::vecS>,
        boost::bidirectionalS,
        DynamicNode,
        DynamicEdge>
// Flagged to kill after the next test 01 jan 2022.
//        boost::property<DynamicNode, boost::vertex_index_t>,
//        boost::property<DynamicEdge, boost::edge_index_t>>
    Graph;


// Flagged to kill after the next test 01 jan 2022.
// Useful types :-)
//property_map<Graph, capacity_t>::type capacity
//        = get(capacity_t(), G);
//property_map<Graph, flow_t>::type flow
//        = get(flow_t(), G);
//
// those can be accessed according to what is explained  here:
// https://www.boost.org/doc/libs/1_77_0/libs/graph/doc/using_adjacency_list.html
//
//typedef boost::property_map<Graph, double DynamicEdge::*>::type DynamicEdgeMap;
//typedef boost::iterator_property_map<std::vector<double>::iterator, DynamicEdgeMap> DynamicEdgeCentralMap;
//typedef boost::iterator_property_map<std::vector<int>::iterator, IndexMap> CentralMap;

typedef boost::property_map<Graph, boost::vertex_index_t>::const_type IndexMap;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
typedef boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
typedef boost::property_map<Graph, boost::edge_owner_t>::const_type EdgeOwnerMap;
typedef boost::property_map<Graph, boost::vertex_local_t>::const_type LocalVertexMap;
typedef boost::property_map<Graph, boost::vertex_global_t>::const_type GlobalVertexMap;

// Flagged to kill after the next test 01 jan 2022.
//using boost::adaptors::transformed;
//#include <boost/graph/adjacency_list.hpp>


class CommonGraphObjectClass{
    /*
     * This  is the general container that all the graphs inherit from :-)
     *
     *                          FUNCTIONS FOR CONFIGURATION
     * Initialization:
     *              the initialization sets the values of the nodes and edges, and it's
     *              expected use is for the setting of the initial values previous to simulations.
     *
     *                          FUNCTIONS FOR DEBUGGING
     * showVertex:
     *              shows the number of nodes and edges each processor 'sees',
     *              and all the nodes are iterated and their central value is reported.
     *              The expected use is confirming the number of required nodes are
     *              those as requested, and also confirming the correct setting of
     *              initial values for simple cases e.g. all nodes equal values
     * showEdges:
     *              shows how many nodes are owned by each processor. This is somehow
     *              complementary to showVertex as in the boost graph distributed framework
     *              processors can 'see' more nodes than those they can modify, thus checking
     *              the sum of nodes with a clear ownership is part of a debugging inspection.
     *              The expected value is that NRequestedNodes = sum(owned nodes)
     * reportNProcs:
     *              shows how many processors are active according to the graph boost interface.
     *              the expected value is that NGraphProcessors = MPI Processors.
     */
    public:
        void showVertex(Graph & g);
        void showEdges(Graph & g);
        void reportNProcs(Graph & g);
        void reportNodes(Graph &g);
        void Initialization(std::vector<std::pair<double, double>> X0_W,
                            double J, Graph & g, unsigned int N);
};


#endif //CPPPROJCT_GENERALGRAPH_H