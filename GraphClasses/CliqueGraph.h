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

    // TODO: move this to CommonGraphObjectClass
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    CliqueGraphObject(unsigned long num_nodes) :
            E((unsigned long)  ((((long) num_nodes) * (long) ((long) num_nodes - (long) 1)) / 2)),
            N(num_nodes),
            g(N) {};

    void build();
};

#endif //CPPPROJCT_CLIQUEGRAPH_H
