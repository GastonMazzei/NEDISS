//
// Created by m4zz31 on 26/10/21.
//
#include  <iostream>
#include <Eigen/Sparse>
#include <string>
#include "../utils/error.h"

#include <omp.h>
#include <vector>
#include <algorithm>

// Parallel & Distributed libraries
#include <boost/mpi.hpp>
#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
//#include <boost/mpi/communicator.hpp>
#include <boost/graph/distributed/named_graph.hpp>
#include <boost/graph/parallel/process_group.hpp>
#include <boost/assert.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/parallel/caching_property_map.hpp>
#include <boost/graph/parallel/algorithm.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/parallel/process_group.hpp>


// Sequential libraries
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/named_graph.hpp>

// Common libraries
#include <boost/graph/graph_traits.hpp>




//--------------------MORE IMPORTS WE ARE NOT USING!--------------------
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/properties.hpp>
//#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/iteration_macros.hpp>
//#include <boost/graph/distributed/concepts.hpp>
//#include <boost/iterator/transform_iterator.hpp>
//#include <boost/property_map/property_map.hpp>
//#include <boost/graph/adjacency_iterator.hpp>
//#include <boost/property_map/parallel/distributed_property_map.hpp>
//#include <boost/property_map/parallel/local_property_map.hpp>
//#include <boost/graph/parallel/detail/property_holders.hpp>
//#include <boost/mpl/if.hpp>
//#include <boost/type_traits/is_same.hpp>
//#include <boost/assert.hpp>
//#include <boost/limits.hpp>
//#include <boost/graph/parallel/properties.hpp>
//#include <boost/graph/parallel/distribution.hpp>
//#include <boost/graph/parallel/algorithm.hpp>
//#include <boost/graph/distributed/selector.hpp>
//#include <boost/graph/parallel/process_group.hpp>
//#include <boost/pending/container_traits.hpp>
//#include <boost/serialization/base_object.hpp>
//#include <boost/mpi/datatype.hpp>
//#include <boost/pending/property_serialize.hpp>
//#include <boost/graph/distributed/unsafe_serialize.hpp>
//// Named vertices
//#include <boost/graph/distributed/shuffled_distribution.hpp>
//-------------------end of more imports---------------


// SMALL WORLD https://www.boost.org/doc/libs/1_74_0/libs/graph/doc/small_world_generator.html
// BGL custom vertex properties https://www.boost.org/doc/libs/1_77_0/libs/graph/doc/using_adjacency_list.html
// Migration to distributed adjacency list https://www.boost.org/doc/libs/1_77_0/libs/graph_parallel/doc/html/distributed_adjacency_list.html
// An exterior property map for the edges should be a double that is used in integration https://www.boost.org/doc/libs/1_55_0/libs/graph/doc/using_property_maps.html
// Example using map and stuff https://github.com/mmccoo/nerd_mmccoo/blob/master/boost_properties/adj_lists.cxx
// Boost general tour: https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/graph/doc/quick_tour.html
// Boost (sequential) adjacency list: https://www.boost.org/doc/libs/1_77_0/libs/graph/doc/using_adjacency_list.html
// Examples of parallel graph usage are (1) https://github.com/boostorg/graph_parallel/tree/master/example
// and (2) https://github.com/boostorg/graph_parallel/tree/master/test

#ifndef CPPPROJCT_GRAPH_H
#define CPPPROJCT_GRAPH_H

using namespace std;
using namespace Eigen;
using namespace boost;


//--------------- SEQUENTIAL ADJACENCY LISTS
//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
// a very humble definition of an edge and a vertex:
// typedef pair<int, int> Edge; typedef graph_traits<Graph>::vertex_descriptor Vertex;



// Dynamic Nodes are the elements we want to evolve over time
// it should have:
// (1) its value:
//      - 1 dimensional for Kuramoto
//      - 3 dimensional for Rossler Oscillators
//      - something still unclear for the d(A)_t \propto \laplace^\alpha (A) model
// (2) NOT IMPLEMENTED: an N-dimensional array with its own properties: characteristic frequency for example
// (3) NOT IMPLEMENTED: a memory capacity that  potentially can remember N_neighbors * msg
struct DynamicNode {
    DynamicNode() { }
    DynamicNode(int value): value(value) { };
    int value;
    // Serialization support is required!
    template<typename Archiver> /*version is const unsigned int*/
    void serialize(Archiver& ar, const unsigned int /*version*/) {
        ar & value;
    }
};






// Edges should have  a (double precission) value which will always
// account for some sort of "interaction"
struct DynamicEdge {
    DynamicEdge() { }
    DynamicEdge(double value) : value(value) { }
    double value;
    // Serialization support is required!
    template<typename Archiver>
    void serialize(Archiver& ar, const unsigned int /*version*/) {
        ar & value;
    }
};




// -----------named nodes were disabled for parallelization with iterators
//
//
//Enabling named vertices
//namespace boost { namespace graph {
//        template<>
//        struct internal_vertex_name<DynamicNode>
//        {
//            typedef multi_index::member<DynamicNode, int, &DynamicNode::name> type;
//        };
//}}
//// Enabling vertex creation in absence of previous existence,
//// in particular avoiding raising an exception
//namespace boost { namespace graph {
//        template<>
//        struct internal_vertex_constructor<DynamicNode>
//        {
//            typedef vertex_from_name<DynamicNode> type;
//        };
//}}



// INTERESTING STRUCTURE FOR A GENERALIZATION
// after https://github.com/mousebird/boost/blob/master/libs/graph_parallel/test/distributed_betweenness_centrality_test.cpp
//#ifdef CSR
//typedef compressed_sparse_row_graph<directedS, no_property, WeightedEdge,
//                                      no_property, distributedS<mpi_process_group> >
//    Graph;
//  typedef compressed_sparse_row_graph<directedS, no_property, WeightedEdge>
//    seqGraph;
//#else
//typedef adjacency_list<vecS,
//        distributedS<mpi_process_group, vecS>,
//        directedS,
//        no_property,
//        property<edge_weight_t, int> > Graph;
//
//typedef adjacency_list<vecS, vecS, directedS, no_property,
//        property<edge_weight_t, int> > seqGraph;
//#endif


typedef adjacency_list<vecS, distributedS<graph::distributed::mpi_process_group, vecS>,
        bidirectionalS, DynamicNode, DynamicEdge> Graph;



class GraphObject {
public:
    Graph g;
    unsigned long N;
    unsigned long E;
    GraphObject();
    GraphObject(int indicated_type, unsigned long num_nodes);
    GraphObject(int indicated_type, unsigned long num_nodes, double probability);

    // a map that could give the value of the nodes but is not working!
    //
    //    typedef property_map<Graph, int DynamicNode::*>::type
//            road_length = get(&DynamicNode::value  , g);
//    one can access the length of any given road based on its edge descriptor e with
//    the expression get(road_length, e),
//    regardless of which process stores the edge

    // Functions for displaying the structure!
    // Potentially drawable if pipelined in graphviz format.
    void showVertex();
    void showEdges();


};



#endif //CPPPROJCT_GRAPH_H


