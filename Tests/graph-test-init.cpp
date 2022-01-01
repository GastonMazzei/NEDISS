//
// Created by m4zz31 on 29/10/21.
//
#include "graph-test-init.h"

// TODO: expand to the remaining graph topologies

void graph_tests_init(int TOPOLOGY, unsigned int SEED, unsigned long N){

    unsigned long K;
    double p;
    if ((TOPOLOGY == 2) || (TOPOLOGY == 3)) {
        p = (double) std::stod(std::getenv("proba"));
    }
    if ((TOPOLOGY == 3)) {
        K = (unsigned long) std::stoul(std::getenv("kneigh"));
    }

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

    } else if (TOPOLOGY == 2) {

        // Erdos Renyi Network
        reproductibility_lock(SEED);
        ErdosRenyiGraphObject G3(N, p);
        test_graph_init<200, ErdosRenyiGraphObject>(G3, "ErdosRenyi");

    } else if (TOPOLOGY == 3) {

        // Small World Network
        reproductibility_lock(SEED);
        SmallWorldGraphObject G4(N, K, p);
        test_graph_init<200, SmallWorldGraphObject>(G4, "SmallWorld");

    }
}