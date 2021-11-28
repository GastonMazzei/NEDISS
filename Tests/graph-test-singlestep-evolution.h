//
// Created by m4zz31 on 31/10/21.
//

#ifndef CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H
#define CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H

#include "test-imports.h"

template <int T, typename GRAPHTYPE, int BATCH, typename DIFFEQ>
void test_graph_singlestep_evolution(GRAPHTYPE &G, 
				                    std::string name,
                                     SolverConfig &SOLVER) {


    // Print in command what test is it
    adsync_message<T>(msg_prev + "'test_" + name + "_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message<T>(msg_prev + "'preparing "+name+" graph for singlestep evolution'", G.g);
    G.build();
    G.Initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier<T>(msg_post + "'preparing ring graph for singlestep evolution'", G.g);


    // Show the total number of created nodes
    adsync_message_barrier<T>(msg_prev + "'reportNodes'", G.g);
    G.reportNodes(G.g);
    adsync_message_barrier<T>(msg_post + "'reportNodes'", G.g);


    // Show nodes
    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    unsigned long NVtot = boost::num_vertices(G.g);
    CommunicationHelper ComHelper(G.g);
    ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
    IntegrationHelper IntHelper(NVtot);
    LayeredSolverHelper LayHelper(NVtot);
    MappingHelper MapHelper(G.g);

    ReferenceContainer REF(ParHelper,
                           ComHelper,
                           G.g,
                           IntHelper,
                           MapHelper,
                           LayHelper,
                           NVtot);


    if (SOLVER.s == 0) {
        GeneralSolver<DIFFEQ, EulerSolver<DIFFEQ>> S_eu("eu",SOLVER.d) ;
        S_eu.SetT0(0);
        S_eu.SetStep(0.01);

        // Test several eq evolutions with Euler
        adsync_message<T>(msg_prev + "'single_eq_evolution' with euler (1 of 3)", G.g);
        single_evolution<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);
        LayHelper.built = true;
        S_eu.EvolveTime();
        adsync_message_barrier<T>(msg_post + "'single_eq_evolution' with euler (1 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_eq_evolution' with euler (2 of 3)", G.g);
        //single_evolution<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);
        single_evolution2<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);

        S_eu.EvolveTime();
        adsync_message_barrier<T>(msg_post + "'single_eq_evolution' with euler (2 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_eq_evolution' with euler (3 of 3)", G.g);
        //single_evolution<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);
        single_evolution2<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);

        S_eu.EvolveTime();
        adsync_message_barrier<T>(msg_post + "'single_eq_evolution' with euler (3 of 3)", G.g);

    } else if (SOLVER.s == 1) {

        GeneralSolver<DIFFEQ, RungeKuttaSolver<DIFFEQ>> S_rk("rk",SOLVER.d, &(SOLVER.P[0])) ;
        S_rk.SetT0(0);
        S_rk.SetStep(0.01);

        // Test several eq evolutions with RungeKutta
        adsync_message<T>(msg_prev + "'single_eq_evolution' with runge kutta (1 of 3)", G.g);
        single_evolution<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, REF, G.N);

        LayHelper.built = true;
        S_rk.EvolveTime();
        adsync_message_barrier<T>(msg_post + "'single_eq_evolution' with runge kutta (1 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_eq_evolution' with runge kutta (2 of 3)", G.g);
        //single_evolution<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, REF, G.N);
        single_evolution2<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, REF, G.N);

        S_rk.EvolveTime();
        adsync_message_barrier<T>(msg_post + "'single_eq_evolution' with runge kutta (2 of 3)", G.g);
        adsync_message<T>(msg_prev + "'single_eq_evolution' with runge kutta (3 of 3)", G.g);
        //single_evolution<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, REF, G.N);
        single_evolution2<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, REF, G.N);

        S_rk.EvolveTime();
        adsync_message_barrier<T>(msg_post + "'single_eq_evolution' with runge kutta (3 of 3)", G.g);
    }

    // Show nodes
    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    // Print in command what test is it and inform that it has ended
    adsync_message<T>(msg_post + "'test_" + name + "_graph_singlestep_evolution'", G.g);
}



template <int BATCH, typename EQCLASS>
void graph_test_singlestep_evolution_helper(unsigned int SEED, unsigned long N, SolverConfig &SOLVER, int TOPOLOGY){

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
        RingGraphObject G(N);
        test_graph_singlestep_evolution<100, RingGraphObject, BATCH, EQCLASS>(G, "Ring", SOLVER);

    } else if (TOPOLOGY == 1) {
        // Clique Network
        CliqueGraphObject G(N);
        test_graph_singlestep_evolution<100, CliqueGraphObject, BATCH, EQCLASS>(G, "Clique", SOLVER);

    } else if (TOPOLOGY == 2) {
        // Erdos Renyi Network
        ErdosRenyiGraphObject G(N, p);
        test_graph_singlestep_evolution<100, ErdosRenyiGraphObject, BATCH, EQCLASS>(G, "ErdosRenyi", SOLVER);

    } else if (TOPOLOGY == 3) {
        // Small World Network
        SmallWorldGraphObject G(N, K, p);
        test_graph_singlestep_evolution<100, SmallWorldGraphObject, BATCH, EQCLASS>(G, "SmallWorld", SOLVER);

    } else error_report("[ERROR] Requested topology does not exist!\n");

}



template <int BATCH>
void graph_tests_singlestep_evolution(unsigned int SEED, unsigned long N, SolverConfig &SOLVER, int TOPOLOGY, int EQNUMBER){

    reproductibility_lock(SEED);
    if (EQNUMBER == 0){
        graph_test_singlestep_evolution_helper<BATCH, NoiselessKuramoto>(SEED, N, SOLVER, TOPOLOGY);
    } else if (EQNUMBER == 1) {
        graph_test_singlestep_evolution_helper<BATCH, LinearTestEquation>(SEED, N, SOLVER, TOPOLOGY);
    } else {
        printf("[FATAL] Required Equation does not exist!\n");
        std::cout<<std::flush;
        exit(1);
    }
}



#endif //CPPPROJCT_GRAPH_TEST_SINGLESTEP_EVOLUTION_H
