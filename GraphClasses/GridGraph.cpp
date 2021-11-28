//
// Created by m4zz31 on 28/11/21.
//

#include "GridGraph.h"
#include "../Utils/adequate_synchronization.h"

using namespace boost;
using namespace std;

void GridGraphObject::build(){
//    if (process_id(g.process_group()) == 0) {
//        int i=0;
//        add_edge(vertex(0, g), vertex(N-1, g), g);
////#pragma omp parallel for
//        for (int i=0; i<N-1; i++){
//            add_edge(vertex(i,g), vertex(i+1,g), g);
//        }
//    }
    printf("Grid constructor not implmented!!!");std::cout<<std::flush;
    exit(1);
    adsync_synchronization_barrier<0>("Grid constructor", g);
}

