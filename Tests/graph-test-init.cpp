//
// Created by m4zz31 on 29/10/21.
//
#include "graph-test-init.h"
#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../GraphClasses/RingGraph.h"


#include "../Utils/reproductibility.h"



void graph_tests_init(unsigned int SEED, int N, double p){
    // ---CONSTRUCTOR AND INITIALIZATION DEBUGGING---

    // Clique Network
    reproductibility_lock(SEED);
    CliqueGraphObject G1(N);
    test_graph_init<200, CliqueGraphObject>(G1, "Clique");

    // Ring Network
    reproductibility_lock(SEED);
    RingGraphObject G2(N);
    test_graph_init<200, RingGraphObject>(G2, "Ring");


    // Erdos Renyi Network
    reproductibility_lock(SEED);
    ErdosRenyiGraphObject G3(N, p);
    test_graph_init<200, ErdosRenyiGraphObject>(G3, "ErdosRenyi");
}