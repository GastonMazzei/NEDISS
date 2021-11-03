//
// Created by m4zz31 on 29/10/21.
//
#include "graph-test-init.h"
#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../GraphClasses/RingGraph.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/global_standard_messages.h"
#include "../Utils/reproductibility.h"

void test_clique_graph_init(int N){
    //                  N
    CliqueGraphObject G(N);

    // Print in command that the test is Clique
    adsync_message(msg_prev + "'test_clique_graph_init'", G.g);

    // Build it
    adsync_message(msg_prev + "'build_clique'", G.g);
    G.build();
    adsync_message_barrier(msg_post + "'build_clique'", G.g);

    // Show the number of created nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Show the number of processors
    adsync_message(msg_prev + "'reportNProcs'", G.g);
    //G.reportNProcs(G.g);
    adsync_message_barrier(msg_post + "'reportNProcs'", G.g);

    // Show the edges
    adsync_message(msg_prev + "'showEdges'", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges'", G.g);

    // Initialize for constant kuramoto
    adsync_message(msg_prev + "'kuramoto_initialization (constant values)'", G.g);
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier(msg_post + "'kuramoto_initialization (constant values)'", G.g);

    // Show the number of created nodes again
    adsync_message(msg_prev + "'showVertex' (post kuramoto constant values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex' (post kuramoto constant values)", G.g);

    // Show the edges again
    adsync_message(msg_prev + "'showEdges' (post kuramoto constant values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges' (post kuramoto constant values)", G.g);

    // Initialize for varied kuramoto
    //--------INITIALIZATION that is not parallelized yet
    std::vector<std::pair<double, double>> X0_W;
    for (int i=0; i<G.N; i++){
            X0_W.push_back({
                std::sin(3.14 * ((double) i) / 13),
                1/((double) G.N) * (double) i
            });
    }
    adsync_message(msg_prev + "'kuramoto_initialization (varied values)'", G.g);
    G.kuramoto_initialization(X0_W, 5.67, G.g, G.N);
    adsync_message_barrier(msg_post + "'kuramoto_initialization (varied values)'", G.g);

    // Show the number of created nodes again
    adsync_message(msg_prev + "'showVertex' (post kuramoto varied values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex' (post kuramoto varied values)", G.g);

    // Show the edges again
    adsync_message(msg_prev + "'showEdges' (post kuramoto varied values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges' (post kuramoto varied values)", G.g);

    // Print in command that the test is Clique
    adsync_message(msg_post + "'test_clique_graph_init'", G.g);
};


void test_erdosRenyi_graph_init(int N, double p){
    //                      N, P
    ErdosRenyiGraphObject G(N, p);
    G.build();

    // Print in command that the test is Erdos Renyi
    adsync_message(msg_prev + "'test_erdosRenyi_graph_init'", G.g);

    // Show the number of created nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Show the number of processors
    adsync_message(msg_prev + "'reportNProcs'", G.g);
    //G.reportNProcs(G.g);
    adsync_message_barrier(msg_post + "'reportNProcs'", G.g);

    // Show the edges
    adsync_message(msg_prev + "'showEdges'", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges'", G.g);

    // Initialize for constant kuramoto
    adsync_message(msg_prev + "'kuramoto_initialization (constant values)'", G.g);
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier(msg_post + "'kuramoto_initialization (constant values)'", G.g);

    // Show the number of created nodes again
    adsync_message(msg_prev + "'showVertex' (post kuramoto constant values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex' (post kuramoto constant values)", G.g);

    // Show the edges again
    adsync_message(msg_prev + "'showEdges' (post kuramoto constant values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges' (post kuramoto constant values)", G.g);

    // Initialize for varied kuramoto
    //--------INITIALIZATION that is not parallelized yet
    std::vector<std::pair<double, double>> X0_W;
    for (int i=0; i<G.N; i++){
        X0_W.push_back({
                               std::sin(3.14 * ((double) i) / 13),
                               1/((double) G.N) * (double) i
                       });
    }
    adsync_message(msg_prev + "'kuramoto_initialization (varied values)'", G.g);
    G.kuramoto_initialization(X0_W, 5.67, G.g, G.N);
    adsync_message_barrier(msg_post + "'kuramoto_initialization (varied values)'", G.g);

    // Show the number of created nodes again
    adsync_message(msg_prev + "'showVertex' (post kuramoto varied values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex' (post kuramoto varied values)", G.g);

    // Show the edges again
    adsync_message(msg_prev + "'showEdges' (post kuramoto varied values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges' (post kuramoto varied values)", G.g);

    // Print in command that the test is Erdos Renyi
    adsync_message(msg_post + "'test_erdosRenyi_graph'", G.g);
};

void test_ring_graph_init(int N){
    //                  N
    RingGraphObject G(N);

    // Print in command that the test is Ring
    adsync_message(msg_prev + "'test_ring_graph_init'", G.g);

    // Build it
    adsync_message(msg_prev + "'build_ring'", G.g);
    G.build();
    adsync_message_barrier(msg_post + "'build_ring'", G.g);

    // Show the nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Show the number of processors
    adsync_message(msg_prev + "'reportNProcs'", G.g);
    //G.reportNProcs(G.g);
    adsync_message_barrier(msg_post + "'reportNProcs'", G.g);

    // Show the edges
    adsync_message(msg_prev + "'showEdges'", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges'", G.g);

    // Initialize for constant kuramoto
    adsync_message(msg_prev + "'kuramoto_initialization (constant values)'", G.g);
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier(msg_post + "'kuramoto_initialization (constant values)'", G.g);

    // Show the number of created nodes again
    adsync_message(msg_prev + "'showVertex' (post kuramoto constant values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex' (post kuramoto constant values)", G.g);

    // Show the edges again
    adsync_message(msg_prev + "'showEdges' (post kuramoto constant values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges' (post kuramoto constant values)", G.g);

    // Initialize for varied kuramoto
    //--------INITIALIZATION that is not parallelized yet
    std::vector<std::pair<double, double>> X0_W;
    for (int i=0; i<G.N; i++){
        X0_W.push_back({
                               std::sin(3.14 * ((double) i) / 13),
                               1/((double) G.N) * (double) i
                       });
    }
    adsync_message(msg_prev + "'kuramoto_initialization (varied values)'", G.g);
    G.kuramoto_initialization(X0_W, 5.67, G.g, G.N);
    adsync_message_barrier(msg_post + "'kuramoto_initialization (varied values)'", G.g);

    // Show the number of created nodes again
    adsync_message(msg_prev + "'showVertex' (post kuramoto varied values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex' (post kuramoto varied values)", G.g);

    // Show the edges again
    adsync_message(msg_prev + "'showEdges' (post kuramoto varied values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier(msg_post + "'showEdges' (post kuramoto varied values)", G.g);

    // Print in command that the test is Ring
    adsync_message(msg_post + "'test_ring_graph_init'", G.g);
};


void graph_tests_init(unsigned int SEED, int N, double p){
    // ---CONSTRUCTOR AND INITIALIZATION DEBUGGING---

    // Clique Network
    //reproductibility_lock(SEED);
    //test_clique_graph_init(N);

    // Ring Network
    //reproductibility_lock(SEED);
    //test_ring_graph_init(N);

    // Erdos Renyi Network
    reproductibility_lock(SEED);
    test_erdosRenyi_graph_init(N, p);
}