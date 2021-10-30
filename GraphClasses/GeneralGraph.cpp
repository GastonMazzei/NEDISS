//
// Created by m4zz31 on 29/10/21.
//

#include "GeneralGraph.h"


using namespace boost;
using namespace std;




void CommonGraphObjectClass::reportNProcs(Graph &g) {
    if (process_id(g.process_group()) == 0) {
        std::cout <<
                  "N processes  is: " <<
                  num_processes(boost::graph::distributed::mpi_process_group()) <<
                  std::endl;
    }
}


void CommonGraphObjectClass::showVertex(Graph &g) {
    // shows all vertex
    auto vpair = vertices(g);
    for(auto v=vpair.first; v!=vpair.second; v++) {
        std::cout << "node w/ value " << g[*v].value << std::endl;
    }
}

void CommonGraphObjectClass::showEdges(Graph &g) {
// show all edges
    auto epair = edges(g);
    for(auto e=epair.first; e!=epair.second; e++) {
        std::cout << "edge w/value " << g[*e].value << std::endl;
    }
}




// ITERATE OVER ALL EDGES LOCALLY IN PARALLEL
// suited for initializing values ;-)
//
//typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap;
//typedef property_map<Graph, vertex_local_t>::const_type LocalMap;
//OwnerMap owner = get(vertex_owner, g);
//LocalMap local = get(vertex_local, g);
//
//unsigned int MY_NUM = process_id(g.process_group());
//unsigned int MY_NODES = num_vertices(g);
//
//cout << " I Have " << MY_NODES << " NODES!" << endl;
//auto vs = vertices(g);
//int counter = 0;
//for (auto vit = vs.first; vit != vs.second; ++vit) {
//counter ++;
//}
//cout << "I have finished travelling over " << counter << " nodes ... " << endl;
//
//      HOW TO INITIALIZE VALUES
//    Graph::vertex_descriptor v = *vertices(g).first;
//    g[v].value = 13;
//    g[v].params.push_back(1);
//    g[v].params.push_back(2);
//    g[v].params.push_back(3);
//    Graph::edge_descriptor e = *out_edges(v, g).first;
//    g[e].value = 3;


// ITERATE OVER ALL THE VERTICES AND ITS NEIGHBORS, unclear if it accesses nonlocal nodes!
//
// Other maybe useful iterators over neighbors are:
// https://www.boost.org/doc/libs/1_52_0/libs/graph/doc/adjacency_iterator.html
//graph_traits<adjacency_list>::out_edge_iterator
//graph_traits<adjacency_list>::adjacency_iterator
//    auto vs = vertices(g);
//    //
//    for (auto vit = vs.first; vit != vs.second; ++vit) {
//        auto neighbors = boost::adjacent_vertices(*vit, g);
//        for (auto nit = neighbors.first; nit != neighbors.second; ++nit)
//            std::cout << *vit << ' ' << *nit << std::endl;
//    }


// ACCESS AND CHANGE NODE AND EDGE VALUES!
//
//      HOW TO INITIALIZE VALUES
//    Graph::vertex_descriptor v = *vertices(g).first;
//    g[v].value = 13;
//    g[v].params.push_back(1);
//    g[v].params.push_back(2);
//    g[v].params.push_back(3);
//    Graph::edge_descriptor e = *out_edges(v, g).first;
//    g[e].value = 3;



// FRAMEWORK THAT IMPLEMENTS MPI AUTOMATIC COMMUNICATIONS SO THAT
// EVERY NODE AND EDGE ARE ACCESSIBLE
//
// ITERATE OVER ALL NODES
//    typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
//    vertex_iterator v, v_end;

// ITERATE OVER ALL EDGES
//typedef graph_traits<Graph>::edge_iterator edge_iterator;
//vertex_iterator e, e_end;

//    for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
//        cout <<
//            "I am process N: " <<
//            process_id(g.process_group()) <<
//            " And this node, which has a local index of: " <<
//            get(local, *v) <<
//            " and is owned by me (Process N: " <<
//            get(owner, *v) <<
//            " ) " <<
//            "Then the graph's value is: " <<
//            g[*v].value ;
//            g[*v].value = 4;
//            cout << g[*v].value <<
//            //"has a 'global' index of: " <<
//            //get(centrality, *v) <<
//            endl;
//    }