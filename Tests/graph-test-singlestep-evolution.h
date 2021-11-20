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

// THIS USES ONLY THE KURAMOTO EQUATION!!!

template <int T, typename GRAPHTYPE, int BATCH>
void test_graph_singlestep_evolution(GRAPHTYPE &G, std::string name,
                                     CommunicationHelper &ComHelper,
                                     ParallelHelper &ParHelper,
                                     IntegrationHelper &IntHelper,
                                     MappingHelper &MapHelper,
                                     int SOLVER) {
    // Print in command what test is it
    adsync_message<T>(msg_prev + "'test_" + name + "_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message<T>(msg_prev + "'preparing "+name+" graph for singlestep evolution'", G.g);
    G.build();
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


    if (SOLVER == 0) {
        GeneralSolver<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>> S_eu("eu",1) ;

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


    } else if (SOLVER == 1) {

        // Please note that:
        // DifferentialEquation has

        double t[4] = {1,2,2,1};
        GeneralSolver<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>> S_rk("rk",4, t) ;

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
    }

    // Show nodes
    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    // Print in command what test is it
    adsync_message<T>(msg_post + "'test_" + name + "_graph_singlestep_evolution'", G.g);
}



template <int BATCH>
void graph_tests_singlestep_evolution(unsigned int SEED, unsigned long N, double p, int SOLVER, int TOPOLOGY){

    reproductibility_lock(SEED);

    if (TOPOLOGY == 0) {
        // Ring Network
        RingGraphObject G(N);

        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        MappingHelper MapHelper(G.g);

        test_graph_singlestep_evolution<100, RingGraphObject, BATCH>(G, "Ring",
                                                                     ComHelper, ParHelper,
                                                                     IntHelper, MapHelper, SOLVER);

    } else if (TOPOLOGY == 1) {
        // Clique Network
        CliqueGraphObject G(N);

        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        MappingHelper MapHelper(G.g);

        test_graph_singlestep_evolution<100, CliqueGraphObject, BATCH>(G, "Clique",
                                                                           ComHelper, ParHelper,
                                                                           IntHelper, MapHelper, SOLVER);

    } else if (TOPOLOGY == 2) {
        // Erdos Renyi Network
        ErdosRenyiGraphObject G(N, p);

        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        MappingHelper MapHelper(G.g);


        test_graph_singlestep_evolution<100, ErdosRenyiGraphObject, BATCH>(G, "ErdosRenyi",
                                                                           ComHelper, ParHelper,
                                                                           IntHelper, MapHelper, SOLVER);
    } else error_report("[ERROR] Requested topology does not exist!\n");
}

#endif //CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H
