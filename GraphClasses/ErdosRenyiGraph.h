//
// Created by m4zz31 on 26/10/21.
//
#ifndef CPPPROJCT_ERDOSRENYIGRAPH_H
#define CPPPROJCT_ERDOSRENYIGRAPH_H

#include "GeneralGraph.h"

class ErdosRenyiGraphObject : public CommonGraphObjectClass {
public:

    unsigned long N;
    unsigned long E;
    Graph g;

    // Unique to this class :-)
    typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
    boost::minstd_rand gen;
    double Proba;

    // TODO: move this to CommonGraphObjectClass
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    ErdosRenyiGraphObject(unsigned long num_nodes, double probability) :
            E(0), // Actually the expected edge number is combinatorial(num_nodes,2) * probability
            N(num_nodes),
            Proba(probability),
            g(ERGen(gen, num_nodes, probability, false), ERGen(), num_nodes) {};

    // TODO: move this to CommonGraphObjectClass
    void build(){}; // empty function for inter-class compatibility :-)
};
#endif //CPPPROJCT_ERDOSRENYIGRAPH_H
