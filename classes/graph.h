//
// Created by m4zz31 on 26/10/21.
//
#include  <iostream>
#include <Eigen/Sparse>
#include <string>
#include "../utils/error.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <boost/mpi.hpp>
#include <boost/graph/use_mpi.hpp>

//#include <boost/graph/distributed/mpi_process_group.hpp>
//#include <boost/graph/erdos_renyi_generator.hpp>

using namespace boost;


// SMALL WORLD https://www.boost.org/doc/libs/1_74_0/libs/graph/doc/small_world_generator.html
// BGL custom vertex properties https://www.boost.org/doc/libs/1_77_0/libs/graph/doc/using_adjacency_list.html
// Migration to distributed adjacency list https://www.boost.org/doc/libs/1_77_0/libs/graph_parallel/doc/html/distributed_adjacency_list.html
// An exterior property map for the edges should be a double that is used in integration https://www.boost.org/doc/libs/1_55_0/libs/graph/doc/using_property_maps.html
// Example using map and stuff https://github.com/mmccoo/nerd_mmccoo/blob/master/boost_properties/adj_lists.cxx


#ifndef CPPPROJCT_GRAPH_H
#define CPPPROJCT_GRAPH_H

using namespace std;
using namespace Eigen;




// more info here: https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/graph/doc/quick_tour.html
class GraphObject {
public:
    typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
    // choosing if list, vec or set: we chose vec
    // https://www.boost.org/doc/libs/1_77_0/libs/graph/doc/using_adjacency_list.html
//    typedef adjacency_list<vecS,
//            distributedS<parallel::mpi::bsp_process_group, vecS>,
//    bidirectionalS> Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef pair<int, int> Edge;
    Graph g;
    GraphObject(int indicated_type, int num_vertices);
    void showVertex();
    void showEdges();
};




#endif //CPPPROJCT_GRAPH_H
//template <class Iter, class ID>
//typename std::iterator_traits<Iter>::value_type
//get(const iterator_map<Iter,ID>& i,
//    typename boost::property_traits<ID>::key_type key)
//{
//    return i.m_iter[i.m_id[key]];
//}
//template <class Iter, class ID>
//void
//put(const iterator_map<Iter,ID>& i,
//    typename boost::property_traits<ID>::key_type key,
//    const typename std::iterator_traits<Iter>::value_type& value)
//{
//    i.m_iter[i.m_id[key]] = value;
//}
//template <class Iter, class ID>
//typename std::iterator_traits<Iter>::reference
//at(const iterator_map<Iter,ID>& i,
//   typename boost::property_traits<ID>::key_type key)
//{
//    return i.m_iter[i.m_id[key]];
//}
