//
// Created by m4zz31 on 29/10/21.
//

#include "GeneralGraph.h"
#include "../Utils/error.h"
#include <cmath>
#include "mpi.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/memory_management.h"
#include "../Utils/msleep.h"
#include <iomanip>

using namespace boost;
using namespace std;

// Please include in your radar the following bug: (potential segfault when one process has 0 nodes)
// https://github.com/boostorg/graph_parallel/issues/18


void CommonGraphObjectClass::reportNProcs(Graph &g) {
    if (process_id(g.process_group()) == 0) {
        std::cout <<
                  "N processes  is: " <<
                  num_processes(boost::graph::distributed::mpi_process_group()) <<
                  std::endl;
    }
    adsync_barrier<0>(); // This barrier (attempts) to fix a bug:
    //                  the program never ends if this function is called w/out barrier
}




void CommonGraphObjectClass::showVertex(Graph &g) {
    unsigned int MY_NUM = process_id(g.process_group());
    unsigned int MY_NODES = num_vertices(g);
    unsigned int MY_EDGES = num_edges(g);

    cout << " I Have " << MY_NODES << " NODES! and " << MY_EDGES << " EDGES!" << endl;
    vertex_iterator v, v_end;
    for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
        if (g[*v].params.size() > 0) {
            std::cout << "node w/ value " << std::setprecision (15)  << g[*v].value << " and first param " << g[*v].params[0] << std::endl;
        } else {
            std::cout << "node w/ value " <<  std::setprecision (15)  << g[*v].value << std::endl;
        }
    }
}

void CommonGraphObjectClass::reportNodes(Graph &g){
    vertex_iterator v, v_end;
    int counter = 0;
    for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
        if (get(get(vertex_owner, g), *v) == process_id(g.process_group())){
            counter++;
        }
    }
    std::cout << "I am Proc N'" << process_id(g.process_group()) <<
              " and, according to 'boost::num_vertices', I report having " <<
              boost::num_vertices(g) <<
              " nodes! After iterating 'boost::vertices' I noted I own " <<
              counter <<
              " nodes." <<
              std::endl;
}




void CommonGraphObjectClass::showEdges(Graph &g) {
    typedef boost::iterator_property_map<std::vector<int>::iterator, IndexMap> CentralMap;
    edge_iterator e, e_end;
    for (boost::tie(e, e_end) = edges(g); e != e_end; ++e) {
        std::cout << "edge w/value " << g[*e].value << std::endl;
        //std::cout << get(&DynamicEdge::value, g) << std::endl;
    }
}


void CommonGraphObjectClass::kuramoto_initialization(std::vector<pair<double, double>> X0_W, double J, Graph & g, unsigned int N){
//    // Whatever Graph Architecture is recieved is populated with a Frequency W_i and Edge weights J

    // 1)
//    // Actually this clause should be changed: we would want to directly recieve W[N] split
//    // into the number of processors, thus wanting W at each processor to be either W[1] or W[N/Nprocs]
    if ((X0_W.size()>1) && (N != X0_W.size())) {
        error_report("[error] Number of frequencies to initialize the Kuramoto model should be either 1 or N");
    } else if (X0_W.size()==1){
        for (int i=0; i<N-1; i++){
            X0_W.push_back(X0_W[0]);
        }
    };


    // 2)
    // Definition of various required types

    OwnerMap owner = get(vertex_owner, g);
    vertex_iterator v, v_end;
    edge_iterator e, e_end;

    // 3)
    // split the data across processors :-)
    static unsigned int NV = num_vertices(g);
    mpi::communicator world;
    unsigned int NPROCS = world.size();
    unsigned int BALANCE = std::floor(((double) N) / ((double) NPROCS));
    unsigned int MY_NUM = process_id(g.process_group());
    pair<double, double> myvals[NV];
    int j = 0;
    for (j=0; j<BALANCE; j++){
        myvals[j] = X0_W[MY_NUM * BALANCE + j];
    }
    if (NV == (BALANCE + 1)){
        myvals[BALANCE] = X0_W[N-1];
    }


    // 4) Iterate along edges and nodes
    j = 0;
    for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
        if (MY_NUM == get(owner, *v)) {
            if (g[*v].params.size()==0) {
                g[*v].params.push_back(myvals[j].second);
            } else {
                g[*v].params[0] = myvals[j].second;
            }
            g[*v].value = myvals[j].first;
            ++j;
        };
    }
    for (boost::tie(e, e_end) = edges(g); e != e_end; ++e) {
        g[*e].value = J;
    };
};








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

