//
// Created by m4zz31 on 29/10/21.
//


#include "GeneralGraph.h"
#include "CliqueGraph.h"
#include "../utils/adequate_synchronization.h"


using namespace std;
using namespace boost;


void CliqueGraphObject::build_clique(){
    typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
    //typedef property_map<Graph, int DynamicNode::*>::type LocalMap;
    //LocalMap local = get(&DynamicNode::value, g);
    //get(g, *v).value = 8;
    if (process_id(g.process_group()) == 0) {
        vertex_iterator v, v_end, v_inner;
        for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
            v_inner = v;
            while (++v_inner != v_end) {
                add_edge(*v, *v_inner, g);
            }
        }
    }
    adsync_synchronization_barrier("Clique constructor", g);
}




