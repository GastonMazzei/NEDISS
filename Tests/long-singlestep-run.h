//
// Created by m4zz31 on 7/11/21.
//

#ifndef CPPPROJCT_LONG_SINGLESTEP_RUN_H
#define CPPPROJCT_LONG_SINGLESTEP_RUN_H


#include "../Utils/adequate_synchronization.h"
#include "../Utils/global_standard_messages.h"
#include "../Solvers/EulerSolver.h"
#include "../Solvers/RungeKuttaSolver.h"
#include "../DifferentialEquations/NoiselessKuramoto.h"
#include "../GraphClasses/GraphFunctions.h"
#include "../Utils/HelperClasses.h"
#include "../Communication/CommunicationFunctions.h"

// THIS USES ONLY THE KURAMOTO EQUATION!!!

template <int T, typename GRAPHTYPE, int BATCH>
void test_long_singlestep_run(GRAPHTYPE &G, std::string name,
                                     CommunicationHelper &ComHelper,
                                     ParallelHelper &ParHelper,
                                     IntegrationHelper &IntHelper,
                                     MappingHelper &MapHelper,
                                     int NRUNS, int SOLVER) {

    // Print in command what test is it
    adsync_message<T>(msg_prev + "'test_" + name + "_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message<T>(msg_prev + "'preparing " + name + " graph for singlestep evolution'", G.g);
    G.build();

    G.kuramoto_initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier<T>(msg_post + "'preparing ring graph for singlestep evolution'", G.g);

    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    if (SOLVER == 0) {
        GeneralSolver<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>> S_eu("eu", 1);
        for (int i = 0; i < NRUNS; i++) {
            // Test several kuramoto evolutions with Euler
            single_evolution<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>, BATCH>(G.g, S_eu, ComHelper, ParHelper,
                                                                                       IntHelper, MapHelper, G.N);
        }
    } else if (SOLVER == 1) {
        double t[4] = {1, 2, 2, 1};
        GeneralSolver<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>> S_rk("rk", 4, t);
        for (int i = 0; i < NRUNS; i++) {
            // Test several kuramoto evolutions with Euler
            single_evolution<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>, BATCH>(G.g, S_rk, ComHelper, ParHelper,
                                                                                       IntHelper, MapHelper, G.N);
        }
    } else error_report("Requested solver does not exist\n");

    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);
}


template <int BATCH>
void central_test_long_singlestep_run(unsigned int SEED, unsigned long N, double p, int NRUNS,  int SOLVER, int TOPOLOGY){
    // Ring Network
    reproductibility_lock(SEED);

    if (TOPOLOGY == 0){
        RingGraphObject G(N);
        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        MappingHelper MapHelper(G.g);
        test_long_singlestep_run<100, RingGraphObject, BATCH>(G, "Ring",
                                                              ComHelper, ParHelper,
                                                              IntHelper, MapHelper,
                                                              NRUNS, SOLVER);
    } else if (TOPOLOGY == 1) {
        CliqueGraphObject G(N);
        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        MappingHelper MapHelper(G.g);
        test_long_singlestep_run<100, CliqueGraphObject, BATCH>(G, "Clique",
                                                                ComHelper, ParHelper,
                                                                IntHelper, MapHelper,
                                                                NRUNS, SOLVER);
    } else if (TOPOLOGY == 2) {
        ErdosRenyiGraphObject G(N, p);
        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        MappingHelper MapHelper(G.g);
        test_long_singlestep_run<100, ErdosRenyiGraphObject, BATCH>(G, "ErdosRenyi",
                                                                    ComHelper, ParHelper,
                                                                    IntHelper, MapHelper,
                                                                    NRUNS, SOLVER);
    } else error_report("Requested ring topology does not exist!\n");

};

#endif //CPPPROJCT_LONG_SINGLESTEP_RUN_H
