//
// Created by m4zz31 on 28/11/21.
//

#ifndef CPPPROJCT_GRIDGRAPH_H
#define CPPPROJCT_GRIDGRAPH_H


#include "GeneralGraph.h"

class GridGraphObject : public CommonGraphObjectClass {
public:
    unsigned long N;
    unsigned long E;
    Graph g;
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    GridGraphObject(unsigned long num_nodes) :
            E(num_nodes), // E =/= 2 * N as it is bidirectional ;-) so it's just N
            N(num_nodes),
            g(N) {};
    void build(); // Ring constructor
};



#endif //CPPPROJCT_GRIDGRAPH_H
