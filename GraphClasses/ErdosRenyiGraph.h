//
// Created by m4zz31 on 26/10/21.
//
#ifndef CPPPROJCT_ERDOSRENYIGRAPH_H
#define CPPPROJCT_ERDOSRENYIGRAPH_H

#include "GeneralGraph.h"

class ErdosRenyiGraphObject : public CommonGraphObjectClass {
public:
    // WE SHOULD CHANGE PROBABILITY FOR NUMBER OF EDGES!
    typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
    boost::minstd_rand gen;
    unsigned long N;
    unsigned long E;
    double Proba;
    Graph g;
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    ErdosRenyiGraphObject(unsigned long num_nodes, double probability) :
            E((unsigned long) (probability * (double) num_nodes)), // THIS IS WRONG!
            N(num_nodes),
            Proba(probability),
            g(ERGen(gen, num_nodes, probability, false), ERGen(), num_nodes) {};

    void build(){}; // empty function for inter-class compatibility :-)
};
#endif //CPPPROJCT_ERDOSRENYIGRAPH_H
