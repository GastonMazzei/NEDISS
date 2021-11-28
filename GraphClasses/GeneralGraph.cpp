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
    }
}



void CommonGraphObjectClass::Initialization(std::vector<pair<double, double>> X0_W, double J, Graph & g, unsigned int N){

    static unsigned int NV = num_vertices(g);
    pair<double, double> myvals[NV];
    OwnerMap owner = get(vertex_owner, g);
    vertex_iterator v, v_end;
    edge_iterator e, e_end;
    unsigned int MY_NUM = process_id(g.process_group());

    if (X0_W.size()>1) {
        PRINTF_DBG("SIZE >1\n");
        if (N == X0_W.size()){
            PRINTF_DBG("SIZE IS N\n");
            // Procedure to select the NNodesProcessor_i with i the current processor
            // for the case in which we are given not 'our' values but those of the entire graph
            mpi::communicator world;
            unsigned int NPROCS = world.size();
            unsigned int BALANCE = (unsigned int) std::floor(((double) N) / ((double) NPROCS));
            bool is_case_below = false;
            bool are_we_below = false;
            if ((N % NPROCS != 0) && (MY_NUM + 1 <= N % NPROCS)) {
                BALANCE++;
                are_we_below = true;
            }
            if (N % NPROCS != 0) {
                is_case_below = true;
            }
            int begin=0,end=0;
            for (int j=0; j<MY_NUM; ++j){
                if (is_case_below){
                    if (are_we_below){
                        begin += (int) BALANCE;
                    } else if (j+1 <= N % NPROCS) {
                        begin += (int) (BALANCE + 1);
                    } else {
                        begin += (int) BALANCE;
                    }
                } else {
                    begin += BALANCE;
                }
            }
            end = begin + BALANCE;
            assert(NV == BALANCE);
            for (int j=begin; j<end; j++){
                myvals[j - begin] = X0_W[j];
            }
        } else if (X0_W.size() != NV) {
            error_report("[error] Length of the params vector to initialize the model should be either 1, NNodesTotal or NNodesProcessor_i");
        } else {
            PRINTF_DBG("SIZE IS NV\n");
            for (int j=0; j<X0_W.size(); j++){
                myvals[j] = X0_W[j];
            }
        }
    } else if (X0_W.size()==1){
        PRINTF_DBG("SIZE IS 1\n");
        for (int j=0; j<NV; j++){
            myvals[j] = X0_W[0];
        }
    } else if (X0_W.size() == 0) {
        error_report("[error] A zero-length parameter vector was passed to initialization.");
    }



    // Iterate along nodes
    int j = 0;
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

    // Iterate along edges
    for (boost::tie(e, e_end) = edges(g); e != e_end; ++e) {
        g[*e].value = J;
    };
};


