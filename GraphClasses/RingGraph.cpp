//
// Created by m4zz31 on 31/10/21.
//
#include "RingGraph.h"
#include "../utils/adequate_synchronization.h"

using namespace boost;
using namespace std;

void RingGraphObject::build_ring(){
    typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
    if (process_id(g.process_group()) == 0) {
        vertex_iterator v, v_end, v_inner;
        tie(v, v_end) = vertices(g);
        add_edge(*v, *(v_end-1), g); // Connect the first one to the last one, i.e. S1 topology
        while (v != (v_end-2)) {
            v_inner = v;
            v_inner ++;
            add_edge(*v, *v_inner, g);
            v += 2;
            add_edge(*v, *v_inner, g);
        }
        add_edge(*v, *(v+1), g); // Connect the last two together,
        //  as applying here the 3-items connectivity algorithm from the loop leads to segmentation fault
        // because there is no g[*(v+2)] ;-)
    }
    adsync_synchronization_barrier("Ring constructor", g);
}

