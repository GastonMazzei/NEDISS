//
// Created by m4zz31 on 29/10/21.
//


#include "GeneralGraph.h"
#include "CliqueGraph.h"
#include "../Utils/adequate_synchronization.h"


using namespace boost;


void CliqueGraphObject::build(){
    if (process_id(g.process_group()) == 0) {
        for (int i = 0; i < N; i++) {
#pragma omp for
            for (int j = i + 1; j < N; j++) {
                add_edge(vertex(i, g), vertex(j, g), g);
            }
        }
    }
    adsync_synchronization_barrier<0>("Clique constructor", g);
}




