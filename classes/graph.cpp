//
// Created by m4zz31 on 26/10/21.
//

#include "graph.h"
#include <boost/random/linear_congruential.hpp>

using namespace std;
using namespace boost;

void GraphObject::showVertex(){
//    auto vpair = vertices(g);
//    for (auto iter=vpair.first; iter!=vpair.second; iter++){
//        cout << "Vertex is: " << *iter << endl;
//    }

}

void GraphObject::showEdges(){
//    auto epair = edges(g);
//    for(auto iter=epair.first; iter!=epair.second; iter++) {
//        std::cout << "edge " << source(*iter, g) << " - " << target(*iter, g) << std::endl;
//    }
}



template<typename RandomGenerator, typename Graph>
class sorted_erdos_renyi_iterator
{
public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::pair<vertices_size_type, vertices_size_type> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef void difference_type;

    sorted_erdos_renyi_iterator();
    sorted_erdos_renyi_iterator(RandomGenerator& gen, vertices_size_type n,
                                double probability = 0.0, bool allow_self_loops = false);

    // Iterator operations
    reference operator*() const;
    pointer operator->() const;
    sorted_erdos_renyi_iterator& operator++();
    sorted_erdos_renyi_iterator operator++(int);
    bool operator==(const sorted_erdos_renyi_iterator& other) const;
    bool operator!=(const sorted_erdos_renyi_iterator& other) const;
};



GraphObject::GraphObject(int indicated_type, unsigned long num_nodes) {
    // indicated_type 0: Undirected & Fully Connected (i.e. Clique)
    // indicated_type 1,2,... others PEND IMPLEMENT
    //
    // Please note that
    // std::minstd_rand is std::linear_congruential_engine<std::uint_fast32_t, 48271, 0, 2147483647>
    //
    // refs of this function are:
    // [1] https://www.boost.org/doc/libs/1_77_0/libs/graph_parallel/doc/html/distributed_adjacency_list.html
    N = num_nodes;
    if (indicated_type != 0) {
        error_report("The indicated type does not exist (line 25 in graph.h)");
    }
    // This is Erdos Renyi, as per the docs [1]
    typedef sorted_erdos_renyi_iterator<minstd_rand, Graph> ERIter;
    minstd_rand gen;
    Graph q(ERIter(gen, N, 0.15, false), ERIter(), N);
    g = q;


    // This is Clique: it should be checked as it is probable that all the processes are lading
    // the same data :-( [1]
//#pragma omp parallel for collapse(2)
//    for (int i = 0; i < num_nodes; i++) {
//        for (int j = 0; j < num_nodes - 1; j++) {
//            if (i != j) {
//                add_edge(1, 2, DynamicEdge(1.), g);
//            };
//        }
//    }
}



