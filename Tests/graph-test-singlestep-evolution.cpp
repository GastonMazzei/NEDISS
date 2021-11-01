//
// Created by m4zz31 on 31/10/21.
//
#include "graph-test-singlestep-evolution.h"
#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../Solvers/GeneralSolvers.h"
#include "../GraphClasses/RingGraph.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/global_standard_messages.h"







void test_clique_graph_singlestep_evolution(int N) {
    //                  N
    CliqueGraphObject G(N);
}

void test_erdosRenyi_graph_singlestep_evolution(int N, double p) {
    //                  N
    ErdosRenyiGraphObject G(N, p);
}

void test_ring_graph_singlestep_evolution(int N) {
    //                  N
    RingGraphObject G(N);
    // Print in command that the test is Ring
    adsync_message(msg_prev + "'test_ring_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message(msg_prev + "'preparing ring graph for singlestep evolution'", G.g);
    G.build_ring();
    Solver s = SolverSelector("eul");
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier(msg_post + "'preparing ring graph for singlestep evolution'", G.g);


    // TEST GOES HERE!
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler", G.g);
    G.single_kuramoto_evolution(G.g, s);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler", G.g);

    // Print in command that the test is Ring
    adsync_message(msg_post + "'test_ring_graph_singlestep_evolution'", G.g);
}