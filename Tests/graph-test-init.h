//
// Created by m4zz31 on 30/10/21.
//

#ifndef CPPPROJCT_GRAPH_TEST_INIT_H
#define CPPPROJCT_GRAPH_TEST_INIT_H

#include <string>
#include <vector>
#include <cmath>

#include "test-imports.h"


#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <boost/graph/graphviz.hpp>

void graph_tests_init(int TOPOLOGY, unsigned int SEED, unsigned long N = 4);

// THIS USES ONLY THE KURAMOTO EQUATION!!!
template <int T, typename GRAPHTYPE>
void test_graph_init(GRAPHTYPE &G, std::string name){

    // Print in command which test is
    adsync_message<T>(msg_prev + "'test_" + name + "_graph_init'", G.g);

    // Build it
    adsync_message<T>(msg_prev + "'build_" + name + "'", G.g);
    G.build();
    adsync_message_barrier<T>(msg_post +  "build_" + name, G.g);

    // Show the total number of created nodes
    adsync_message_barrier<T>(msg_prev + "'reportNodes'", G.g);
    G.reportNodes(G.g);
    adsync_message_barrier<T>(msg_post + "'reportNodes'", G.g);

    // Show the created nodes
    adsync_message_barrier<T>(msg_prev + "'showVertex'", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex'", G.g);

    // Show the edges
    adsync_message_barrier<T>(msg_prev + "'showEdges'", G.g);
    G.showEdges(G.g);
    adsync_message_barrier<T>(msg_post + "'showEdges'", G.g);

    // Initialize for constant vals
    adsync_message_barrier<T>(msg_prev + "'Initialization' (constant values)'", G.g);
    G.Initialization({{12.345, 6.78}}, 3.14, G.g, G.N);
    adsync_message_barrier<T>(msg_post + "'Initialization' (constant values)'", G.g);

    // Show the number of created nodes again
    adsync_message_barrier<T>(msg_prev + "'showVertex' (post eq constant values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex' (post eq constant values)", G.g);

    // Show the edges again
    adsync_message_barrier<T>(msg_prev + "'showEdges' (post eq constant values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier<T>(msg_post + "'showEdges' (post eq constant values)", G.g);

    // Initialize for varied kuramoto
    //--------INITIALIZATION that is not parallelized yet
    std::vector<std::pair<double, double>> X0_W;
    for (int i=0; i<G.N; i++){
        X0_W.push_back({
                               std::abs(std::sin(3.14 * ((double) i) / 13)),
                               1/((double) G.N) * (double) i
                       });
    }
    adsync_message_barrier<T>(msg_prev + "'Initialization' (varied values)'", G.g);
    G.Initialization(X0_W, 5.67, G.g, G.N);
    adsync_message_barrier<T>(msg_post + "'Initialization' (varied values)'", G.g);

    // Show the number of created nodes again
    adsync_message_barrier<T>(msg_prev + "'showVertex' (post eq varied values)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex' (post eq varied values)", G.g);

    // Show the number of created nodes again
    adsync_message_barrier<T>(msg_prev + "'showVertex' (post eq varied values)(again, to confirm the vertex iterator conserves the order)", G.g);
    G.showVertex(G.g);
    adsync_message_barrier<T>(msg_post + "'showVertex' (post eq varied values)(again, to confirm the vertex iterator conserves the order)", G.g);


    // Show the edges again
    adsync_message_barrier<T>(msg_prev + "'showEdges' (post eq varied values)", G.g);
    G.showEdges(G.g);
    adsync_message_barrier<T>(msg_post + "'showEdges' (post eq varied values)", G.g);

    // Print in command what test is :-)
    adsync_message<T>(msg_post + "'test_" + name + "_graph_init'", G.g);

    const char* dot = "graphviz_test_new.dot";
    boost::write_graphviz(dot, G.g, boost::make_label_writer(get(&DynamicNode::value, G.g)));

//    template <class Name>
//    class label_writer {
//    public:
//        label_writer(Name _name) : name(_name) {}
//        template <class VertexOrEdge>
//        void operator()(std::ostream& out, const VertexOrEdge& v) const {
//            out << "[label=\"" << name[v] << "\"]";
//        }
//    private:
//        Name name;
//    };
//    A function to conveniently create this writer is provided:
//
//    template < class Name >
//    label_writer<Name>
    //auto boost::make_label_writer(Name n);


};



#endif //CPPPROJCT_GRAPH_TEST_INIT_H
