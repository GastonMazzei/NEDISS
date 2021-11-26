//
// Created by m4zz31 on 26/11/21.
//

#ifndef CPPPROJCT_SMALLWORLDGRAPH_H
#define CPPPROJCT_SMALLWORLDGRAPH_H

#include "GeneralGraph.h"
#include <boost/graph/small_world_generator.hpp>

class SmallWorldGraphObject : public CommonGraphObjectClass {
public:
    typedef boost::small_world_iterator<boost::minstd_rand, Graph> SWGen;
    boost::minstd_rand gen;
    unsigned long N;
    unsigned long K;
    unsigned long E;
    double Proba;
    Graph g;
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    SmallWorldGraphObject(unsigned long num_nodes, unsigned long kNeighbors, double probability) :
            E(kNeighbors * num_nodes),
            K(kNeighbors),
            N(num_nodes),
            Proba(probability),
            g(SWGen(gen, num_nodes, kNeighbors, probability, false), SWGen(), num_nodes) {};

    void build(){}; // empty function for inter-class compatibility :-)
};


#endif //CPPPROJCT_SMALLWORLDGRAPH_H
