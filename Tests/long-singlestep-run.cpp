//
// Created by m4zz31 on 7/11/21.
//
#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"
#include "../GraphClasses/RingGraph.h"
#include "../Utils/reproductibility.h"
#include "long-singlestep-run.h"
#include "../Communication/CommunicationFunctions.h"
#include "../GraphClasses/GraphFunctions.h"


void central_test_long_singlestep_run(unsigned int SEED, int N, double p, int NRUNS){
    // Ring Network
    reproductibility_lock(SEED);
    RingGraphObject G2(N);

    // helpers instantiated here just temporaly :-)
    unsigned long NVtot = boost::num_vertices(G2.g);
    CommunicationHelper ComHelper(G2.g);
    ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
    IntegrationHelper IntHelper(NVtot);
    MappingHelper MapHelper(G2.g);

    test_long_singlestep_run<100, RingGraphObject>(G2, "Ring",
                                                      ComHelper, ParHelper,
                                                      IntHelper, MapHelper,
                                                        NRUNS);
};