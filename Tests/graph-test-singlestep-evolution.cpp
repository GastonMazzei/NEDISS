//
// Created by m4zz31 on 31/10/21.
//
#include "graph-test-singlestep-evolution.h"

#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../GraphClasses/RingGraph.h"
#include "../Utils/reproductibility.h"


#include "../Communication/CommunicationFunctions.h"
#include "../GraphClasses/GraphFunctions.h"

void graph_tests_singlestep_evolution(unsigned int SEED, int N, double p){
    // ---SINGLE TIMESTEP EVOLUTION DEBUGGING---

    // Clique Network
//    reproductibility_lock(SEED);
//    CliqueGraphObject G1(N);
//    test_graph_singlestep_evolution<100, CliqueGraphObject>(G1, "Clique");

    // Ring Network
    reproductibility_lock(SEED);
    RingGraphObject G2(N);

    // helpers instantiated here just temporaly :-)
    unsigned long NVtot = boost::num_vertices(G2.g);
    CommunicationHelper ComHelper(G2.g);
    ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
    IntegrationHelper IntHelper(NVtot);
    MappingHelper MapHelper(G2.g);

    test_graph_singlestep_evolution<100, RingGraphObject>(G2, "Ring",
                                                          ComHelper, ParHelper,
                                                          IntHelper, MapHelper);
    //CommunicationHelper H(G2.g);
    //GetMsgFromSomeone(H,G2.g);




//    // Erdos Renyi Network
//    reproductibility_lock(SEED);
//    ErdosRenyiGraphObject G3(N, p);
//    test_graph_singlestep_evolution<100, ErdosRenyiGraphObject>(G3, "ErdosRenyi");
}