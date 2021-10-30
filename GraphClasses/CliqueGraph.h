//
// Created by m4zz31 on 29/10/21.
//

#ifndef CPPPROJCT_CLIQUEGRAPH_H
#define CPPPROJCT_CLIQUEGRAPH_H


#include "GeneralGraph.h"


class CliqueGraphObject : public CommonGraphObjectClass {
public:
    unsigned long N;
    unsigned long E;
    Graph g;
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    CliqueGraphObject(unsigned long num_nodes) :
            E((num_nodes * (num_nodes - 1)) / 2),
            N(num_nodes),
            g(N) {};

    void build_clique();
};

#endif //CPPPROJCT_CLIQUEGRAPH_H
