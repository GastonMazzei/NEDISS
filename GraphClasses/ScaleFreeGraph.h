//
// Created by m4zz31 on 28/11/21.
//

#ifndef CPPPROJCT_SCALEFREEGRAPH_H
#define CPPPROJCT_SCALEFREEGRAPH_H

#include "GeneralGraph.h"
#include <boost/graph/plod_generator.hpp>


// TODO: test this class
class ScaleFreeGraphObject : public CommonGraphObjectClass {
public:
    unsigned long N;
    unsigned long E;
    Graph g;

    // Uniques to this Graph :-)
    typedef boost::plod_iterator<boost::minstd_rand, Graph> SFGen;
    boost::minstd_rand gen;
    double alpha,beta;

    // TODO: move this to CommonGraphObjectClass
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    ScaleFreeGraphObject(unsigned long num_nodes, double A, double B) :
            alpha(A),
            beta(B),
            N(num_nodes),
            g(SFGen(gen, num_nodes, A, B, false), SFGen(), num_nodes) {};

    // TODO: move this to CommonGraphObjectClass
    void build(){};
};


#endif //CPPPROJCT_SCALEFREEGRAPH_H
