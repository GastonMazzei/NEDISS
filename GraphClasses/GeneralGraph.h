//
// Created by m4zz31 on 29/10/21.
//

#ifndef CPPPROJCT_GENERALGRAPH_H
#define CPPPROJCT_GENERALGRAPH_H

#include  <iostream>
#include <Eigen/Sparse>
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


// Dynamic Nodes are the elements we want to evolve over time
// it should have:
// (1) its value:
//      - 1 dimensional for Kuramoto
//      - 3 dimensional for Rossler Oscillators
//      - something still unclear for the d(A)_t \propto \laplace^\alpha (A) model
// (2) params:
//      - an N-dimensional array with its own properties: characteristic frequency for example
// (3) temoral register:
//      - store the result of operations aand keep it for some time
struct DynamicNode {
    double value = 0;
    double temporal_register = 0;
    DynamicNode() = default;
    DynamicNode(double i): value(i) {};
    std::vector<double> params;
    template<typename Archiver> /*version is const unsigned int*/
    void serialize(Archiver& ar, const unsigned int /*version*/) {
        ar & value & params ;
    }
};


// Edges should have  a (double precission) value which will always
// account for some sort of "interaction"
struct DynamicEdge {
    DynamicEdge() = default;
    DynamicEdge(double i): value(i) {};
    double value = 1;
    // Serialization support is required!
    template<typename Archiver>
    void serialize(Archiver& ar, const unsigned int /*version*/) {
        ar & value;
    }
};

// Unclear if it is necessary
typedef DynamicNode DynamicNode;
typedef DynamicEdge DynamicEdge;


// A central object in this work: the "Graph" type
typedef boost::adjacency_list<boost::vecS,
        boost::distributedS<boost::graph::distributed::mpi_process_group, boost::vecS>,
        boost::bidirectionalS,
        DynamicNode,
        DynamicEdge>
//        boost::property<DynamicNode, boost::vertex_index_t>,
//        boost::property<DynamicEdge, boost::edge_index_t>>
    Graph;



// Useful types :-)
//property_map<Graph, capacity_t>::type capacity
//        = get(capacity_t(), G);
//property_map<Graph, flow_t>::type flow
//        = get(flow_t(), G);
//
// those can be accessed according to what is explained  here:
// https://www.boost.org/doc/libs/1_77_0/libs/graph/doc/using_adjacency_list.html
//
typedef boost::property_map<Graph, boost::vertex_index_t>::const_type IndexMap;
//typedef boost::property_map<Graph, double DynamicEdge::*>::type DynamicEdgeMap;
//typedef boost::iterator_property_map<std::vector<double>::iterator, DynamicEdgeMap> DynamicEdgeCentralMap;
typedef boost::iterator_property_map<std::vector<int>::iterator, IndexMap> CentralMap;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
typedef boost::property_map<Graph, boost::vertex_owner_t>::const_type OwnerMap;
typedef boost::property_map<Graph, boost::edge_owner_t>::const_type EdgeOwnerMap;
typedef boost::property_map<Graph, boost::vertex_local_t>::const_type LocalVertexMap;
typedef boost::property_map<Graph, boost::vertex_global_t>::const_type GlobalVertexMap;

// ---------------------------------------------------------------------------------
// AN ANSWER IN STACKOVERFLOW (https://stackoverflow.com/questions/68936738/iterate-over-bundled-properties-in-boost-graph-in-a-direct-way)
// TO THE PROBLEM OF NOT BEING ABLE TO BUILD A VERTEX LIST ADAPTOR AFTER https://www.boost.org/doc/libs/1_77_0/libs/graph_parallel/doc/html/vertex_list_adaptor.html
// ---------------------------------------------------------------------------------
//#include <boost/graph/adjacency_list.hpp>
#include <boost/range/adaptors.hpp>
//using boost::adaptors::transformed;


class CommonGraphObjectClass{
    public:
        //using vertex = typename graph::vertex_descriptor
        void showVertex(Graph & g);
        void showEdges(Graph & g);
        void reportNProcs(Graph & g);
        void reportNodes(Graph &g);
        void Initialization(std::vector<std::pair<double, double>> X0_W, double J, Graph & g, unsigned int N);
};


#endif //CPPPROJCT_GENERALGRAPH_H



//    int world_rank;
//    int world_size;
//
//    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
