//
// Created by m4zz31 on 31/10/21.
//

#ifndef CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H
#define CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H


#include "../Utils/adequate_synchronization.h"
#include "../Utils/global_standard_messages.h"
#include "../Solvers/EulerSolver.h"
#include "../Solvers/RungeKuttaSolver.h"
#include "../DifferentialEquations/NoiselessKuramoto.h"
#include "../GraphClasses/GraphFunctions.h"
#include "../Utils/HelperClasses.h"
#include "../Communication/CommunicationFunctions.h"


#include "graph-test-singlestep-evolution.h"

#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../GraphClasses/RingGraph.h"
#include "../Utils/reproductibility.h"


#include "../Communication/CommunicationFunctions.h"
#include "../GraphClasses/GraphFunctions.h"


template <int T, typename GRAPHTYPE, int BATCH>
void test_graph_singlestep_evolution(GRAPHTYPE &G, std::string name,
                                     CommunicationHelper &ComHelper,
                                     ParallelHelper &ParHelper,
                                     IntegrationHelper &IntHelper,
                                     MappingHelper &MapHelper) {
    // Print in command what test is it
    adsync_message<T>(msg_prev + "'test_" + name + "_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message<T>(msg_prev + "'preparing "+name+" graph for singlestep evolution'", G.g);
    G.build();
    GeneralSolver<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>> S_eu("eu",1) ;
    double t[4] = {1,2,2,1};
    GeneralSolver<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>> S_rk("rk",4, t) ;
    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier<T>(msg_post + "'preparing ring graph for singlestep evolution'", G.g);


    // Show the total number of created nodes
    adsync_message_barrier<T>(msg_prev + "'reportNodes'", G.g);
    G.reportNodes(G.g);
    adsync_message_barrier<T>(msg_post + "'reportNodes'", G.g);


//     Show nodes
    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    if (true) {
        // Test several kuramoto evolutions with Euler
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
        single_evolution<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>, BATCH>(G.g, S_eu,ComHelper, ParHelper,
                                                                            IntHelper, MapHelper, G.N);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
        single_evolution<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>, BATCH>(G.g, S_eu,ComHelper, ParHelper,
                                                                            IntHelper, MapHelper, G.N);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);
        single_evolution<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>, BATCH>(G.g, S_eu,ComHelper, ParHelper,
                                                                            IntHelper, MapHelper, G.N);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);


        // Test several kuramoto evolutions with RungeKutta
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
        single_evolution<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>, BATCH>(G.g, S_rk,ComHelper, ParHelper,
                                                                                 IntHelper, MapHelper, G.N);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
        single_evolution<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>, BATCH>(G.g, S_rk,ComHelper, ParHelper,
                                                                                 IntHelper, MapHelper, G.N);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);
        single_evolution<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>, BATCH>(G.g, S_rk,ComHelper, ParHelper,
                                                                                 IntHelper, MapHelper, G.N);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);
    };

    // Show nodes
    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    // Print in command what test is it
    adsync_message<T>(msg_post + "'test_" + name + "_graph_singlestep_evolution'", G.g);
}



template <int BATCH>
void graph_tests_singlestep_evolution(unsigned int SEED, int N, double p){
    // ---SINGLE TIMESTEP EVOLUTION DEBUGGING---

    // Clique Network
//    reproductibility_lock(SEED);
//    CliqueGraphObject G1(N);
//    test_graph_singlestep_evolution<100, ErdosRenyiGraphObject, BATCH>(G1, "Clique",
//                                                                 ComHelper, ParHelper,
//                                                                 IntHelper, MapHelper);

    // Ring Network
    reproductibility_lock(SEED);
    RingGraphObject G2(N);

    // helpers instantiated here just temporaly :-)
    unsigned long NVtot = boost::num_vertices(G2.g);
    CommunicationHelper ComHelper(G2.g);
    ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
    IntegrationHelper IntHelper(NVtot);
    MappingHelper MapHelper(G2.g);

    test_graph_singlestep_evolution<100, RingGraphObject, BATCH>(G2, "Ring",
                                                          ComHelper, ParHelper,
                                                          IntHelper, MapHelper);




//    // Erdos Renyi Network
//    reproductibility_lock(SEED);
//    ErdosRenyiGraphObject G3(N, p);
//    test_graph_singlestep_evolution<100, ErdosRenyiGraphObject, BATCH>(G3, "ErdosRenyi",
//                                                                 ComHelper, ParHelper,
//                                                                 IntHelper, MapHelper);
}

#endif //CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H
