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

