//
// Created by m4zz31 on 31/10/21.
//
#include "graph-test-singlestep-evolution.h"
#include "../Solvers/EulerSolver.h"
#include "../Solvers/RungeKuttaSolver.h"
#include "../DifferentialEquations/NoiselessKuramoto.h"
#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../GraphClasses/GraphFunctions.h"
#include "../Solvers/GeneralSolver.h"
#include "../GraphClasses/RingGraph.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/global_standard_messages.h"
#include "../Utils/reproductibility.h"





void test_clique_graph_singlestep_evolution(int N) {
    //                  N
    CliqueGraphObject G(N);
    // Print in command that the test is Ring
    adsync_message<0>(msg_prev + "'test_clique_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message(msg_prev + "'preparing clique graph for singlestep evolution'", G.g);
    G.build();
    GeneralSolver<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>> S_eu("eu",1) ;
    double t[4] = {1,2,2,1};
    GeneralSolver<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>> S_rk("rk",4, t) ;
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier(msg_post + "'preparing clique graph for singlestep evolution'", G.g);

    // Show nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Test several kuramoto evolutions with Euler
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);


    // Test several kuramoto evolutions with RungeKutta
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);

    // Show nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Print in command that the test is Ring
    adsync_message(msg_post + "'test_clique_graph_singlestep_evolution'", G.g);

}

void test_erdosRenyi_graph_singlestep_evolution(int N, double p) {
    //                  N
    ErdosRenyiGraphObject G(N, p);
    G.build();

    // Print in command that the test is Ring
    adsync_message(msg_prev + "'test_erdosRenyi_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message(msg_prev + "'preparing erdosRenyi graph for singlestep evolution'", G.g);
    GeneralSolver<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>> S_eu("eu",1) ;
    double t[4] = {1,2,2,1};
    GeneralSolver<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>> S_rk("rk",4, t) ;
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier(msg_post + "'preparing erdosRenyi graph for singlestep evolution'", G.g);

    // Show nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Test several kuramoto evolutions with Euler
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);


    // Test several kuramoto evolutions with RungeKutta
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);

    // Show nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Print in command that the test is Ring
    adsync_message(msg_post + "'test_erdosRenyi_graph_singlestep_evolution'", G.g);

}


void test_ring_graph_singlestep_evolution(int N) {
    //                  N
    RingGraphObject G(N);
    // Print in command that the test is Ring
    adsync_message(msg_prev + "'test_ring_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message(msg_prev + "'preparing ring graph for singlestep evolution'", G.g);
    G.build();
    GeneralSolver<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>> S_eu("eu",1) ;
    double t[4] = {1,2,2,1};
    GeneralSolver<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>> S_rk("rk",4, t) ;
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier(msg_post + "'preparing ring graph for singlestep evolution'", G.g);

    // Show nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Test several kuramoto evolutions with Euler
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);
    single_evolution<NoiselessKuramoto,EulerSolver<NoiselessKuramoto>>(G.g, S_eu);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);


    // Test several kuramoto evolutions with RungeKutta
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
    adsync_message(msg_prev + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);
    single_evolution<NoiselessKuramoto,RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk);
    adsync_message_barrier(msg_post + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);


    // Show nodes
    adsync_message(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier(msg_post + "'showVertex'", G.g);

    // Print in command that the test is Ring
    adsync_message(msg_post + "'test_ring_graph_singlestep_evolution'", G.g);
}


void graph_tests_singlestep_evolution(unsigned int SEED, int N, double p){
    // ---SINGLE TIMESTEP EVOLUTION DEBUGGING---

    // Clique Network
    reproductibility_lock(SEED);
    test_clique_graph_singlestep_evolution(N);

    // Ring Network
    //reproductibility_lock(SEED);
    //test_ring_graph_singlestep_evolution(N);

    // Erdos Renyi Network
    //reproductibility_lock(SEED);
    //test_erdosRenyi_graph_singlestep_evolution(N, p);
}