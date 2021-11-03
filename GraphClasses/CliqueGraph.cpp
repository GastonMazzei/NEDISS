//
// Created by m4zz31 on 29/10/21.
//


#include "GeneralGraph.h"
#include "CliqueGraph.h"
#include "../Utils/adequate_synchronization.h"


using namespace std;
using namespace boost;


void CliqueGraphObject::build(){
    typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
    //typedef property_map<Graph, int DynamicNode::*>::type LocalMap;
    //LocalMap local = get(&DynamicNode::value, g);
    //get(g, *v).value = 8;
    int i=0, j=0;
    if (process_id(g.process_group()) == 0) {
#pragma omp parallel for
        for (int i=0; i<N; i++){
            for (int j=i+1; j<N; j++){
                add_edge(vertex(i,g), vertex(j,g), g);
            }
        }
    }
    adsync_synchronization_barrier<0>("Clique constructor", g);
}




