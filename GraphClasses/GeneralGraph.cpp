//
// Created by m4zz31 on 29/10/21.
//

#include "GeneralGraph.h"



void CommonGraphObjectClass::reportNProcs(Graph &g){
    std::cout << "N processes  is: " << num_processes(boost::graph::distributed::mpi_process_group()) << std::endl;
}


void CommonGraphObjectClass::showVertex(Graph &g){

//    typedef property_map<Graph, vertex_index_t>::const_type IndexMap;
    // instantiate and ?set properties?
    //IndexMap im(get(vertex_index, g));
    //vector_property_map<DynamicNode, IndexMap> props(num_vertices(g), im);
    //using Traits = graph_traits<Graph>;
    //property_map<Graph, DynamicNode>::const_type
    //get(DynamicNode, g, X x);
//    typedef iterator_property_map<std::vector<int>::iterator, IndexMap>
//            CentralityMap;
//    std::vector<int> centralityS(num_vertices(*g), 0);
//    CentralityMap centrality(centralityS.begin(), get(vertex_index, *g));
//    typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
//    vertex_iterator v, v_end;
//    typedef property_map<Graph, vertex_owner_t>::const_type OwnerMap;
//    typedef property_map<Graph, vertex_local_t>::const_type LocalMap;
//    OwnerMap owner = get(vertex_owner, *g);
//    LocalMap local = get(vertex_local, *g);
//    cout << "before loop" << endl;
//    for (boost::tie(v, v_end) = vertices(*g); v != v_end; ++v) {
//        //cout << get(centrality, *v) << endl;
//        cout << "Lap" << endl;
//        //cout << get(owner, *v) << endl;
//        //out << get(local, *v) << endl;
//    }
//    cout << "after loop " << endl;
    //
//    graph_traits<Graph>::vertex_descriptor start = vertex(0, g);
//    property_map<Graph, DynamicNode>::type value =
//            get(vertex_distance, g);
    std::cout << "About to call!" << std::endl;
    for (int i=0; i<num_processes(boost::graph::distributed::mpi_process_group()); i++){
        if (process_id(g.process_group()) == i) {
            // show all your vertex please :-)
            std::cout << "Hey! :-) I am processor " << process_id(g.process_group()) << std::endl;
            std::cout << "The number of nodes I have are: " << num_vertices(g) << std::endl;
        }
    }
}

void CommonGraphObjectClass::showEdges(Graph & g){
//    auto epair = edges(g);
//    for(auto iter=epair.first; iter!=epair.second; iter++) {
//        std::cout << "edge " << source(*iter, g) << " - " << target(*iter, g) << std::endl;
//    }
}