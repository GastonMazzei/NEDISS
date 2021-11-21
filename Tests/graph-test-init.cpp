//
// Created by m4zz31 on 29/10/21.
//
#include "graph-test-init.h"
#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../GraphClasses/RingGraph.h"


#include "../Utils/reproductibility.h"



void graph_tests_init(int TOPOLOGY, unsigned int SEED, unsigned long N, double p){
    // ---CONSTRUCTOR AND INITIALIZATION DEBUGGING---

    if (TOPOLOGY == 0) {

        // Ring Network
        reproductibility_lock(SEED);
        RingGraphObject G2(N);
        test_graph_init<200, RingGraphObject>(G2, "Ring");

    } else if (TOPOLOGY == 1) {

        // Clique Network
        reproductibility_lock(SEED);
        CliqueGraphObject G1(N);
        test_graph_init<200, CliqueGraphObject>(G1, "Clique");

    } else if (TOPOLOGY == 3) {

        // Erdos Renyi Network
        reproductibility_lock(SEED);
        ErdosRenyiGraphObject G3(N, p);
        test_graph_init<200, ErdosRenyiGraphObject>(G3, "ErdosRenyi");

    }
}