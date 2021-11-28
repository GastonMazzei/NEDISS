//
// Created by m4zz31 on 28/11/21.
//

#ifndef CPPPROJCT_SCALEFREEGRAPH_H
#define CPPPROJCT_SCALEFREEGRAPH_H

#include "GeneralGraph.h"
#include <boost/graph/plod_generator.hpp>

class ScaleFreeGraphObject : public CommonGraphObjectClass {
public:
    typedef boost::plod_iterator<boost::minstd_rand, Graph> SFGen;
    boost::minstd_rand gen;
    unsigned long N;
    double alpha,beta;
    unsigned long E;
    Graph g;
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    ScaleFreeGraphObject(unsigned long num_nodes, double A, double B) :
            alpha(A),
            beta(B),
            N(num_nodes),
            g(SFGen(gen, num_nodes, A, B, false), SFGen(), num_nodes) {};

    void build(){}; // empty function for inter-class compatibility :-)
};


#endif //CPPPROJCT_SCALEFREEGRAPH_H
