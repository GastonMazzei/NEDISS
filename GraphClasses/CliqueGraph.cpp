//
// Created by m4zz31 on 29/10/21.
//


#include "GeneralGraph.h"
#include "CliqueGraph.h"



using namespace std;
using namespace boost;


void CliqueGraphObject::build_clique(){
    if (process_id(g.process_group()) == 0) {
//#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            add_vertex(DynamicNode(0), g);
            //put(IndexMap, g);
            cout << "Added one vertex..." << endl;
        }
//#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                add_edge(vertex(i, g), vertex(j, g), g);
            }
        }
    }
    synchronize(g.process_group());
}