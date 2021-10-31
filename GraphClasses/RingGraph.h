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
            E(2 * num_nodes),
            N(num_nodes),
            g(N) {};

    void build_ring();
};







#endif //CPPPROJCT_RINGGRAPH_H
