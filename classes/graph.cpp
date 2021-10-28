//
// Created by m4zz31 on 26/10/21.
//

#include "graph.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>

using namespace std;
using namespace boost;


typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<Graph>::edge_descriptor edge_descriptor;







// STUFF FOR PARALLEL BOT GRAPH
//
//    typedef boost::graph::parallel::process_group_type<Graph>::type
//            process_group_type;
//    process_group_type pg = process_group(g);
//
//    process_group_type::process_id_type id = process_id(pg);
//    process_group_type::process_size_type p = num_processes(pg);





void GraphObject::showVertex(){
//    auto vpair = vertices(g);
//    for (auto iter=vpair.first; iter!=vpair.second; iter++){
//        cout << "Vertex is: " << *iter << endl;
//    }

    typedef property_map<Graph, vertex_index_t>::const_type IndexMap;
    typedef iterator_property_map<std::vector<int>::iterator, IndexMap>
            CentralityMap;

    std::vector<int> centralityS(num_vertices(g), 0);
    CentralityMap centrality(centralityS.begin(), get(vertex_index, g));

    typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
    vertex_iterator v, v_end;

    // We're cheating when we map vertices in g to vertices in sg
    // since we know that g is using the block distribution here
//    typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap;
//    typedef property_map<Graph, vertex_local_t>::const_type LocalMap;
//    OwnerMap owner = get(vertex_owner, g);
//    LocalMap local = get(vertex_local, g);


//    for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
//        cout << get(centrality, *v) << endl;
//        cout << get(owner, *v) << endl;
//        cout << get(local, *v) << endl;
//    }
//
//    graph_traits<Graph>::vertex_descriptor start = vertex(0, g);
//    property_map<Graph, DynamicNode>::type value =
//            get(vertex_distance, g);
}

void GraphObject::showEdges(){
//    auto epair = edges(g);
//    for(auto iter=epair.first; iter!=epair.second; iter++) {
//        std::cout << "edge " << source(*iter, g) << " - " << target(*iter, g) << std::endl;
//    }
}



GraphObject::GraphObject(int indicated_type, unsigned long num_nodes, double probability){
    // ALLOWED: Erdos Renyi (1)
    // NOT IMPLEMENTED YET: ...
    // DISALLOWED: Clique (0)

    if (indicated_type == 0) {
        error_report("Dont initialize Clique using a probability: exiting because user might be confused");
    }

    // if this is Erdos Renyi then
    N = num_nodes;
    E = (long) (probability * (double) N);

    typedef sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
    minstd_rand gen;
#pragma omp parallel
    // I guess this should only change something if each processor implements threading internally
    {
        Graph g(ERGen(gen, N, E, false), ERGen(), N);
    }

    //gen.seed(SEED_INT); reset seed?

    // Note: edges and vertex still lack initialization for their internal values :-)
}

GraphObject::GraphObject(int indicated_type, unsigned long num_nodes) {
    // ALLOWED: Clique (0)
    // NOT IMPLEMENTED YET: ...
    // DISALLOWED: Erdos Renyi (1)

    if (indicated_type == 1) {
        error_report("Erdos Renyi requires a probability");
    }

    // if this is Clique then
    N = num_nodes;
    E = N * (N-1);

    // DOING THIS IN PARALLEL REQUIRES:
    // https://www.boost.org/doc/libs/1_58_0/libs/graph_parallel/doc/html/index.html
    // (accepted answer in) https://stackoverflow.com/questions/30135470/random-access-of-vertices-using-boostgraph


    Graph g(N);
    if (process_id(g.process_group()) == 0) {
//#pragma omp parallel for
            for (int i = 0; i < num_nodes; i++) {
                add_vertex(DynamicNode(0), g);
            }
//#pragma omp parallel for
            for (int i = 0; i < num_nodes; i++) {
                for (int j = i; j < num_nodes; j++) {
                    add_edge(vertex(i, g), vertex(j, g), g);
                }
            }
    }
    synchronize(g.process_group());
}



