//
// Created by m4zz31 on 7/11/21.
//

#ifndef CPPPROJCT_LONG_SINGLESTEP_RUN_H
#define CPPPROJCT_LONG_SINGLESTEP_RUN_H

#include "test-imports.h"


template <int T, typename GRAPHTYPE, int BATCH, typename DIFFEQ>
void test_long_singlestep_run(GRAPHTYPE &G, std::string name,
                                     int NRUNS,
                                    SolverConfig &SOLVER) {

    // Print in command what test is it
    adsync_message<T>(msg_prev + "'test_" + name + "_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message<T>(msg_prev + "'preparing " + name + " graph for singlestep evolution'", G.g);
    G.build();
    double J = 3.14;
    double W0 = 6.78;
    if (std::stoi(std::getenv("EQNUMBER")) == 0){
        // If it is kuramoto, get the coupling strength from the environment variables
        J = (double) std::stod(std::getenv("J"));
    }
    G.Initialization({{12.345, W0}}, J, G.g, G.N);
    adsync_message_barrier<T>(msg_post + "'preparing ring graph for singlestep evolution'", G.g);

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
        GeneralSolver<DIFFEQ, EulerSolver<DIFFEQ>> S_eu("eu", SOLVER.d);
        S_eu.SetT0(0);
        S_eu.SetStep(0.01);
        for (int i = 0; i < NRUNS; i++) {
            if (i==0){
                single_evolution<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);
            } else {
                single_evolution2<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);
            }
            if (i==0)  LayHelper.built = true;
            S_eu.EvolveTime();
        }
    } else if (SOLVER.s == 1) {
        GeneralSolver<DIFFEQ, RungeKuttaSolver<DIFFEQ>> S_rk("rk", SOLVER.d, &SOLVER.P[0]);
        S_rk.SetT0(0);
        S_rk.SetStep(0.01);
        for (int i = 0; i < NRUNS; i++) {
            if (i==0) {
                single_evolution<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, REF, G.N);
            } else {
                single_evolution2<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, REF, G.N);
            }
            if (i==0)  LayHelper.built = true;
            S_rk.EvolveTime();
        }
    } else error_report("Requested solver does not exist\n");

    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);
}


template <int BATCH, typename EQCLASS>
void central_test_long_singlestep_run_helper(unsigned int SEED, unsigned long N, int NRUNS,  SolverConfig &SOLVER, int TOPOLOGY){

    unsigned long K;
    double p;
    if ((TOPOLOGY == 2) || (TOPOLOGY == 3)) {
        p = (double) std::stod(std::getenv("proba"));
    }
    if ((TOPOLOGY == 3)) {
        K = (unsigned long) std::stoul(std::getenv("kneigh"));
    }

    if (TOPOLOGY == 0){
        RingGraphObject G(N);
        test_long_singlestep_run<100, RingGraphObject, BATCH, EQCLASS>(G, "Ring",
                                                                       NRUNS,
                                                                       SOLVER);
    } else if (TOPOLOGY == 1) {
        CliqueGraphObject G(N);
        test_long_singlestep_run<100, CliqueGraphObject, BATCH, EQCLASS>(G, "Clique",
                                                                         NRUNS,
                                                                         SOLVER);
    } else if (TOPOLOGY == 2) {
        ErdosRenyiGraphObject G(N, p);
        test_long_singlestep_run<100, ErdosRenyiGraphObject, BATCH, EQCLASS>(G, "ErdosRenyi",
                                                                             NRUNS,
                                                                             SOLVER);
    } else if (TOPOLOGY == 3) {
        SmallWorldGraphObject G(N, K, p);
        test_long_singlestep_run<100, SmallWorldGraphObject, BATCH, EQCLASS>(G, "SmallWorld",
                                                                             NRUNS,
                                                                             SOLVER);
    } else error_report("Requested ring topology does not exist!\n");

}


template <int BATCH>
void central_test_long_singlestep_run(unsigned int SEED, unsigned long N, int NRUNS,  SolverConfig &SOLVER, int TOPOLOGY, int EQNUMBER){
    // Ring Network
    reproductibility_lock(SEED);

    if (EQNUMBER == 0){
        central_test_long_singlestep_run_helper<BATCH, NoiselessKuramoto>(SEED, N,NRUNS,  SOLVER, TOPOLOGY);
    } else if (EQNUMBER == 1) {
        central_test_long_singlestep_run_helper<BATCH, LinearTestEquation>(SEED, N, NRUNS,  SOLVER, TOPOLOGY);
    } else {
        printf("[FATAL] Required Equation does not exist!\n");
        std::cout<<std::flush;
        exit(1);
    }

};

#endif //CPPPROJCT_LONG_SINGLESTEP_RUN_H
