//
// Created by m4zz31 on 31/10/21.
//
#include "RingGraph.h"
#include "../Utils/adequate_synchronization.h"

using namespace boost;
using namespace std;

void RingGraphObject::build(){
    if (process_id(g.process_group()) == 0) {
        int i=0;
        add_edge(vertex(0, g), vertex(N-1, g), g);
#pragma omp parallel for
        for (int i=0; i<N-1; i++){
            add_edge(vertex(i,g), vertex(i+1,g), g);
        }
    }
    adsync_synchronization_barrier<0>("Ring constructor", g);
}

