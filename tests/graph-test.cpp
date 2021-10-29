//
// Created by m4zz31 on 29/10/21.
//

#include "../GraphClasses/ErdosRenyiGraph.h"
#include "../GraphClasses/CliqueGraph.h"

void test_clique_graph(int N){
    //                  N
    CliqueGraphObject G(N);
    G.showVertex(G.g);
    G.reportNProcs(G.g);
    G.showEdges(G.g);
};


void test_erdosRenyi_graph(int N, double p){
    //                      N, P
    ErdosRenyiGraphObject G(N, p);
    G.showVertex(G.g);
    G.reportNProcs(G.g);
    G.showEdges(G.g);
};