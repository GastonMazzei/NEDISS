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
#include "../DifferentialEquations/LinearTestEquation.h"
#include "../GraphClasses/GraphFunctions.h"
#include "../Utils/HelperClasses.h"
#include "../Communication/CommunicationFunctions.h"

// THIS USES ONLY THE KURAMOTO EQUATION!!!

template <int T, typename GRAPHTYPE, int BATCH, typename DIFFEQ>
void test_long_singlestep_run(GRAPHTYPE &G, std::string name,
                                     CommunicationHelper &ComHelper,
                                     ParallelHelper &ParHelper,
                                     IntegrationHelper &IntHelper,
                                     MappingHelper &MapHelper,
                                     int NRUNS,
                                    LayeredSolverHelper &LayHelper,
                                    SolverConfig &SOLVER) {

    // Print in command what test is it
    adsync_message<T>(msg_prev + "'test_" + name + "_graph_singlestep_evolution'", G.g);

    // Preprocessing
    adsync_message<T>(msg_prev + "'preparing " + name + " graph for singlestep evolution'", G.g);
    G.build();

    G.Initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier<T>(msg_post + "'preparing ring graph for singlestep evolution'", G.g);

    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    if (SOLVER.s == 0) {
        GeneralSolver<DIFFEQ, EulerSolver<DIFFEQ>> S_eu("eu", SOLVER.d);
        S_eu.SetT0(0);
        S_eu.SetStep(0.01);
        for (int i = 0; i < NRUNS; i++) {
            // Test several kuramoto evolutions with Euler
            single_evolution<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, ComHelper, ParHelper,
                                                                   IntHelper, MapHelper, LayHelper, G.N);
            if (i==0)  LayHelper.built = true;
            S_eu.EvolveTime();
        }
    } else if (SOLVER.s == 1) {
        GeneralSolver<DIFFEQ, RungeKuttaSolver<DIFFEQ>> S_rk("rk", SOLVER.d, &SOLVER.P[0]);
        S_rk.SetT0(0);
        S_rk.SetStep(0.01);
        for (int i = 0; i < NRUNS; i++) {
            // Test several kuramoto evolutions with Euler
            single_evolution<DIFFEQ, RungeKuttaSolver<DIFFEQ>, BATCH>(G.g, S_rk, ComHelper, ParHelper,
                                                                       IntHelper, MapHelper, LayHelper, G.N);
            if (i==0)  LayHelper.built = true;
            S_rk.EvolveTime();
        }
    } else error_report("Requested solver does not exist\n");

    adsync_message<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);
}


template <int BATCH, typename EQCLASS>
void central_test_long_singlestep_run_helper(unsigned int SEED, unsigned long N, double p, int NRUNS,  SolverConfig &SOLVER, int TOPOLOGY){

    if (TOPOLOGY == 0){
        RingGraphObject G(N);
        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        LayeredSolverHelper LayHelper(NVtot);
        MappingHelper MapHelper(G.g);
        test_long_singlestep_run<100, RingGraphObject, BATCH, EQCLASS>(G, "Ring",
                                                                       ComHelper, ParHelper,
                                                                       IntHelper, MapHelper,
                                                                       NRUNS, LayHelper, SOLVER);
    } else if (TOPOLOGY == 1) {
        CliqueGraphObject G(N);
        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        LayeredSolverHelper LayHelper(NVtot);
        MappingHelper MapHelper(G.g);
        test_long_singlestep_run<100, CliqueGraphObject, BATCH, EQCLASS>(G, "Clique",
                                                                         ComHelper, ParHelper,
                                                                         IntHelper, MapHelper,
                                                                         NRUNS, LayHelper, SOLVER);
    } else if (TOPOLOGY == 2) {
        ErdosRenyiGraphObject G(N, p);
        // helpers instantiated here just temporaly :-)
        unsigned long NVtot = boost::num_vertices(G.g);
        CommunicationHelper ComHelper(G.g);
        ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
        IntegrationHelper IntHelper(NVtot);
        LayeredSolverHelper LayHelper(NVtot);
        MappingHelper MapHelper(G.g);
        test_long_singlestep_run<100, ErdosRenyiGraphObject, BATCH, EQCLASS>(G, "ErdosRenyi",
                                                                             ComHelper, ParHelper,
                                                                             IntHelper, MapHelper,
                                                                             NRUNS, LayHelper, SOLVER);
    } else error_report("Requested ring topology does not exist!\n");

}


template <int BATCH>
void central_test_long_singlestep_run(unsigned int SEED, unsigned long N, double p, int NRUNS,  SolverConfig &SOLVER, int TOPOLOGY, int EQNUMBER){
    // Ring Network
    reproductibility_lock(SEED);

    if (EQNUMBER == 0){
        central_test_long_singlestep_run_helper<BATCH, NoiselessKuramoto>(SEED, N, p, NRUNS,  SOLVER, TOPOLOGY);
    } else if (EQNUMBER == 1) {
        central_test_long_singlestep_run_helper<BATCH, LinearTestEquation>(SEED, N, p, NRUNS,  SOLVER, TOPOLOGY);
    } else {
        printf("[FATAL] Required Equation does not exist!\n");
        std::cout<<std::flush;
        exit(1);
    }

};

#endif //CPPPROJCT_LONG_SINGLESTEP_RUN_H
