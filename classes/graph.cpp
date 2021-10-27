//
// Created by m4zz31 on 26/10/21.
//
#include <Eigen/Sparse>
#include <string>
#include "graph.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <algorithm>
#include <algorithm>

#include <boost/graph/erdos_renyi_generator.hpp>

using namespace std;
using namespace boost;

void GraphObject::showVertex(){
    auto vpair = vertices(g);
    for (auto iter=vpair.first; iter!=vpair.second; iter++){
        cout << "Vertex is: " << *iter << endl;
    }
}

void GraphObject::showEdges(){
    auto epair = edges(g);
    for(auto iter=epair.first; iter!=epair.second; iter++) {
        std::cout << "edge " << source(*iter, g) << " - " << target(*iter, g) << std::endl;
    }
}

GraphObject::GraphObject(int indicated_type, int num_vertices) {
    // indicated_type 0: Undirected & Fully Connected (i.e. Clique)
    // indicated_type 1,2,... others PEND IMPLEMENT
    //
    // num_vertices... an integer
    //
    if (indicated_type != 0) {
        error_report("The indicated type does not exist (line 25 in graph.h)");
    }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices - 1; j++) {
            if (i != j) {
                add_edge(i, j, g);
            };
        }
    }
}



