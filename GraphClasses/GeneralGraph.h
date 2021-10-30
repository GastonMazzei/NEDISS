//
// Created by m4zz31 on 29/10/21.
//

#ifndef CPPPROJCT_GENERALGRAPH_H
#define CPPPROJCT_GENERALGRAPH_H

#include  <iostream>
#include <Eigen/Sparse>
#include <string>
#include "../utils/error.h"

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
// (3) NOT IMPLEMENTED: a memory capacity that  potentially can remember N_neighbors * msg
struct DynamicNode {
    DynamicNode() = default;
    DynamicNode(int i): value(i) {};
    int value = 0;
    std::vector<int> params;
    template<typename Archiver> /*version is const unsigned int*/
    void serialize(Archiver& ar, const unsigned int /*version*/) {
        ar & value & params ;
    }
};


// Edges should have  a (double precission) value which will always
// account for some sort of "interaction"
struct DynamicEdge {
    DynamicEdge() = default;
    DynamicEdge(int i): value(i) {};
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
        DynamicEdge, boost::vertex_index_t>
    Graph;


// Printing and other stuff
class CommonGraphObjectClass{
    public:
        void showVertex(Graph & g);
        void showEdges(Graph & g);
        void reportNProcs(Graph & g);
};


#endif //CPPPROJCT_GENERALGRAPH_H
