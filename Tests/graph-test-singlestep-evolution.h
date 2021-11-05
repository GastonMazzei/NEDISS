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

template <int T, typename GRAPHTYPE>
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

    // Show nodes
    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    if (true) {
        // Test several kuramoto evolutions with Euler
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
        single_evolution<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>>(G.g, S_eu,ComHelper, ParHelper,
                                                                            IntHelper, MapHelper);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with euler (1 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
        single_evolution<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>>(G.g, S_eu,ComHelper, ParHelper,
                                                                            IntHelper, MapHelper);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with euler (2 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);
        single_evolution<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>>(G.g, S_eu,ComHelper, ParHelper,
                                                                            IntHelper, MapHelper);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with euler (3 of 3)", G.g);


        // Test several kuramoto evolutions with RungeKutta
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
        single_evolution<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk,ComHelper, ParHelper,
                                                                                 IntHelper, MapHelper);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with runge kutta (1 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
        single_evolution<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk,ComHelper, ParHelper,
                                                                                 IntHelper, MapHelper);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with runge kutta (2 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);
        single_evolution<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>>(G.g, S_rk,ComHelper, ParHelper,
                                                                                 IntHelper, MapHelper);
        adsync_message_barrier<T>(msg_post + "'single_kuramoto_evolution' with runge kutta (3 of 3)", G.g);
    };

    // Show nodes
    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    // Print in command what test is it
    adsync_message<T>(msg_post + "'test_" + name + "_graph_singlestep_evolution'", G.g);
}


void graph_tests_singlestep_evolution(unsigned int SEED, int N = 4, double p = 0.4);

#endif //CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H
