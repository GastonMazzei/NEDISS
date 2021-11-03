//
// Created by m4zz31 on 31/10/21.
//

#ifndef CPPPROJCT_RINGGRAPH_H
#define CPPPROJCT_RINGGRAPH_H

#include "GeneralGraph.h"



class RingGraphObject : public CommonGraphObjectClass {
public:
    unsigned long N;
    unsigned long E;
    Graph g;
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    RingGraphObject(unsigned long num_nodes) :
            E(num_nodes), // E =/= 2 * N as it is bidirectional ;-) so it's just N
            N(num_nodes),
            g(N) {};

    void build(); // Ring constructor
//    using vertex = typename Graph::vertex_descriptor;
//    auto bundles() {
//        auto accessor = [map = get(boost::vertex_bundle, g)](vertex v) -> auto& {
//            return map[v];
//        };
//        return vertices(g) | transformed(accessor);
//    }
//    auto bundles() const {
//        auto accessor = [map = get(boost::vertex_bundle, g)](vertex v) -> auto& {
//            return map[v];
//        };
//        return vertices(g) | transformed(accessor);
//    }
};







#endif //CPPPROJCT_RINGGRAPH_H
