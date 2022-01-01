//
// Created by m4zz31 on 26/11/21.
//

#ifndef CPPPROJCT_SMALLWORLDGRAPH_H
#define CPPPROJCT_SMALLWORLDGRAPH_H

#include "GeneralGraph.h"
#include <boost/graph/small_world_generator.hpp>

class SmallWorldGraphObject : public CommonGraphObjectClass {
public:
    /*
     *  The small world graph can be constructed using the Watts-Strogatz
     *  mechanism which imposes the constraint of no self-loops AND no double edges.
     *  The Wikipedia article (https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model)
     *  mentions the constraint N >> K >> log(N), which might have to do with
     *  the mechanism effectively recovering the definition of small world, i.e. that the
     *  typical distance = L propto log(N) in presence of the two previously mentioned constraints
     *
     *  TLDR: if something weird is happening with your results, please notify developers AND
     *  consider the bound K >> log(N).
     *
     *  Quantitative enforcement:
     *  using the base 2 for the log, and defining '>>' as two orders of magnitude, i.e. a factor 4,
     *  then N >> K >> log(N) implies that
     *          1) N >= 16 * log(N), which is enforced with ~ N >= 100
     *          2) a formula for K(N) such that K is equidistant in exponential 2-steps, i.e.
     *             log2(K/N) ~= log2(log2(N)/K) so a safe generator for K is sqrt(log2(N)*N).
     *             If N is increased enough then just for example 100 * log2(N) will suffice
     *             mini-table:
     *             N     K=sqrt(log2(N)*N)  log2(N)   K~100 * log2(N)
     *            1e3        100             10           X
     *            1e4        350             15          1500
     *            1e5       1300             17          2000
     *            1e6       4500             20          2000
     *            1e7      15000             23          2500
     *            1e8      50000             27          2500
     */
    unsigned long N;
    unsigned long E;
    Graph g;

    // Unique to this Graph :-)
    typedef boost::small_world_iterator<boost::minstd_rand, Graph> SWGen;
    boost::minstd_rand gen;
    unsigned long K;
    double Proba;

    // TODO: move this to CommonGraphObjectClass
    unsigned long Procs = num_processes(boost::graph::distributed::mpi_process_group());

    SmallWorldGraphObject(unsigned long num_nodes, unsigned long kNeighbors, double probability) :
            E(kNeighbors * num_nodes),
            K(kNeighbors),
            N(num_nodes),
            Proba(probability),
            g(SWGen(gen, num_nodes, kNeighbors, probability, false), SWGen(), num_nodes) {};

    // TODO: move this to CommonGraphObjectClass
    void build(){};
};


#endif //CPPPROJCT_SMALLWORLDGRAPH_H
