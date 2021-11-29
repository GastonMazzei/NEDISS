//
// Created by m4zz31 on 26/11/21.
//

#ifndef CPPPROJCT_SIMULATIONEXAMPLE_H
#define CPPPROJCT_SIMULATIONEXAMPLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <boost/graph/graphviz.hpp>

#include "../Tests/test-imports.h"

template <int T, typename GRAPHTYPE, int BATCH, typename DIFFEQ>
void run_several_times(GRAPHTYPE &G, std::string name,
//                              CommunicationHelper &ComHelper,
//                              ParallelHelper &ParHelper,
//                              IntegrationHelper &IntHelper,
//                              MappingHelper &MapHelper,
                              int NRUNS,
//                              LayeredSolverHelper &LayHelper,
                              SolverConfig &SOLVER) {


    // Build in case the chosen topology requires building
    G.build();

    // Initialize values
    std::vector<std::pair<double, double>> X0_W;
    double J = 3.14;
    double WMIN=1,WMAX=2;
    if (std::stoi(std::getenv("EQNUMBER")) == 0){
        // If it is kuramoto, get the coupling strength from the environment variables
        J = (double) std::stod(std::getenv("J"));
        WMIN = (double) std::stod(std::getenv("WMIN"));
        WMAX = (double) std::stod(std::getenv("WMAX"));
    }
    std::uniform_real_distribution<double> initialPhase(0, M_PI);
    std::uniform_real_distribution<double> naturalFrequency(WMIN, WMAX);
    std::default_random_engine re;
    for (int i=0; i<G.N; i++){
        double new_phase = initialPhase(re);
        double new_freq = naturalFrequency(re);
        X0_W.push_back({
                        new_phase,new_freq
// deterministic alternative:
//                       std::abs(std::sin(3.14 * ((double) i) / 13)),
//                       1/((double) G.N) * (double) i
                       });
    }

    G.Initialization(X0_W, J, G.g, G.N);

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


    int SAMPLING_FREQ = (int) std::stoi(std::getenv("SAMPLING_FREQ"));
    // Start the simulation!
    if (SOLVER.s == 0) {
        int fileCounter = 0;
        int WRITING_FREQ = SAMPLING_FREQ;
        GeneralSolver<DIFFEQ, EulerSolver<DIFFEQ>> S_eu("eu", SOLVER.d);
        S_eu.SetT0(0);
        S_eu.SetStep(0.01);
        for (int i = 0; i < NRUNS; i++) {
            if (i==0) {
                single_evolution<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);
            } else {
                single_evolution2<DIFFEQ, EulerSolver<DIFFEQ>, BATCH>(G.g, S_eu, REF, G.N);
            }
            if (i==0)  LayHelper.built = true;
            S_eu.EvolveTime();
            if ((i % WRITING_FREQ) == 0) {
                std::string dots = "graphic/program-output/graphviz." + std::to_string(fileCounter) + ".dot";
                char * dot = &dots[0];
                boost::write_graphviz(dot, G.g, boost::make_label_writer(get(&DynamicNode::value, G.g)));
                fileCounter++;
                std::cout << "Writing graphviz into file named" << dots <<std::endl;
            }
            // Evolve the graph? if it were dynamic. Not our case :-).
        }
    } else if (SOLVER.s == 1) {
        int fileCounter = 0;
        int WRITING_FREQ = SAMPLING_FREQ;
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
            if ((i % WRITING_FREQ) == 0) {
                std::string dots = "graphic/program-output/graphviz." + std::to_string(fileCounter) + ".dot";
                char * dot = &dots[0];
                boost::write_graphviz(dot, G.g, boost::make_label_writer(get(&DynamicNode::value, G.g)));
                fileCounter++;
                std::cout << "Writing graphviz into file named" << dots <<std::endl;
            }
            // Evolve the graph? if it were dynamic. Not our case :-).
        }
    } else error_report("Requested solver does not exist\n");

    // End
    printf("\n\n\n[INFO] Simulation has ended successfully!!!\n\n\n");

}


template <int BATCH, typename EQCLASS>
void simulator_helper(unsigned int SEED, unsigned long N, int NRUNS,  SolverConfig &SOLVER, int TOPOLOGY){

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
        run_several_times<100, RingGraphObject, BATCH, EQCLASS>(G, "Ring",
//                                                                       ComHelper, ParHelper,
//                                                                       IntHelper, MapHelper,
                                                                       NRUNS,
//                                                                       LayHelper,
                                                                       SOLVER);
    } else if (TOPOLOGY == 1) {
        CliqueGraphObject G(N);
        run_several_times<100, CliqueGraphObject, BATCH, EQCLASS>(G, "Clique",
//                                                                         ComHelper, ParHelper,
//                                                                         IntHelper, MapHelper,
                                                                         NRUNS,
//                                                                         LayHelper,
                                                                         SOLVER);
    } else if (TOPOLOGY == 2) {
        ErdosRenyiGraphObject G(N, p);
        run_several_times<100, ErdosRenyiGraphObject, BATCH, EQCLASS>(G, "ErdosRenyi",
//                                                                             ComHelper, ParHelper,
//                                                                             IntHelper, MapHelper,
                                                                             NRUNS,
//                                                                             LayHelper,
                                                                             SOLVER);
    } else if (TOPOLOGY == 3) {
        SmallWorldGraphObject G(N, K, p);
        run_several_times<100, SmallWorldGraphObject, BATCH, EQCLASS>(G, "SmallWorld",
//                                                                             ComHelper, ParHelper,
//                                                                             IntHelper, MapHelper,
                                                                             NRUNS,
//                                                                             LayHelper,
                                                                             SOLVER);
    } else error_report("Requested ring topology does not exist!\n");

}


template <int BATCH>
void simulator(unsigned int SEED, unsigned long N, int NRUNS,  SolverConfig &SOLVER, int TOPOLOGY, int EQNUMBER){

    reproductibility_lock(SEED);

    if (EQNUMBER == 0){
        simulator_helper<BATCH, NoiselessKuramoto>(SEED, N,NRUNS,  SOLVER, TOPOLOGY);
    } else if (EQNUMBER == 1) {
        simulator_helper<BATCH, LinearTestEquation>(SEED, N, NRUNS,  SOLVER, TOPOLOGY);
    } else {
        printf("[FATAL] Required Equation does not exist!\n");
        std::cout<<std::flush;
        exit(1);
    }

};


#endif //CPPPROJCT_SIMULATIONEXAMPLE_H
