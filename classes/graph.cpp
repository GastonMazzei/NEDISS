//
// Created by m4zz31 on 26/10/21.
//

#include "graph.h"

using namespace std;
using namespace boost;

void GraphObject::showVertex(){
//    auto vpair = vertices(g);
//    for (auto iter=vpair.first; iter!=vpair.second; iter++){
//        cout << "Vertex is: " << *iter << endl;
//    }

}

void GraphObject::showEdges(){
//    auto epair = edges(g);
//    for(auto iter=epair.first; iter!=epair.second; iter++) {
//        std::cout << "edge " << source(*iter, g) << " - " << target(*iter, g) << std::endl;
//    }
}



GraphObject::GraphObject(int indicated_type, int num_nodes) {
    // indicated_type 0: Undirected & Fully Connected (i.e. Clique)
    // indicated_type 1,2,... others PEND IMPLEMENT
    //
    // num_vertices... an integer
    //
    if (indicated_type != 0) {
        error_report("The indicated type does not exist (line 25 in graph.h)");
    }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < num_nodes - 1; j++) {
            if (i != j) {
                add_edge(1, 2, DynamicEdge(1.), g);
            };
        }
    }
}



