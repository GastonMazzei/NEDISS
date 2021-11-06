//
// Created by m4zz31 on 3/11/21.
//

#include "GraphFunctions.h"


void register_to_value(Graph &g){
    auto vs = vertices(g);
    for (auto v = vs.first; v != vs.second; ++v) {
        g[*v].value = g[*v].temporal_register;
    }
};

void contribute_to_integration(ReferenceContainer &REF){
    // ******************HOW THE INTEGRATION CONTINUES  ;-)*****************
    //
    // Perform the evolution and store the result in the central
    // node's temporal register
//        g[*v].temporal_register = solver.evolve(centralValue,
//                                             centralParams,
//                                             neighborValues,
//                                             edgeValues);
//        clear_vectors(centralParams, neighborValues, edgeValues);
}
