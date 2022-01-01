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

/*
                                     BUG REPORT:
         1)
          asking for N nodes and P (mpi) processors produces a sigfault if P>N
          https://github.com/boostorg/graph_parallel/issues/18
*/



void CommonGraphObjectClass::reportNProcs(Graph &g) {
    /*
     * Report how many processors are active according to the boost graph interface
     */
    if (process_id(g.process_group()) == 0) {
        std::cout <<
                  "N processes  is: " <<
                  num_processes(boost::graph::distributed::mpi_process_group()) <<
                  std::endl;
    }
    adsync_barrier<0>();
}




void CommonGraphObjectClass::showVertex(Graph &g) {
    /*
    * Report how many nodes and edges does each processor can see.
    * Also, the nodes are iterated and their value is shown.
    * It could be useful for checking initialization in simple cases where
    * the idea is to e.g. start being all the same
    */
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
    /*
     * Report how many nodes does each processor own
     */
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
    /*
     * Function to initialize graphs
     *
     *                     NODE INITIALIZATION SUPPORT (X0_W -> Nodes Values)
     *
     *      1) X0_W has length == Total Number of Nodes. The mapping is done using
     *         the vector's indexes as the number of node index, which is assigned
     *         ~randomly by the Boost Graph interface.
     *              comment:
     *                      this should be enough for the study of attractors
     *                      starting in random configurations, but can't produce
     *                      a specific initial value configuration
     *
     *      2) X0_W has length == Total Number of locally owned nodes. This is as
     *         the previous case but specially appropiate if an (either random or
     *         constant) initializer of length total number of nodes cannot fit
     *         comfortably in memory.
     *              comment:
     *                      same as (1).
     *
     *      3) X0_W has length == 1. This is interpreted as a command for initializing
     *         all nodes to the same value.
     *              comment:
     *                      as a constant initializer interface is more comfortable
     *                      than producing a 'number of locally owned nodes' vector
     *                      before calling this function, but this function will build
     *                      it anyway so it shouldn't be confused with a space optimization.
     *
     *      Further Comment:
     *              For initial value conditions problems it is more appropriate to
     *              initialize both the graph itself and the values from a file :-)
     *
     *                     NODE INITIALIZATION SUPPORT (X0_W -> Nodes Values)
     *
     *      1) J is a constant value, and all the edges will have the same value.
     *
     *
     *      Further Comment:
     *              As before, to set specific conditions (e.g. distance-dependent-coupling)
     *              the correct path would be to initialize both the graph and the values
     *              from a file :-)
     *
     *      TODO:
     *            1) Build the graph building + initialization from file
     *            2) Refactor this function as to accept either vectors for random J's,
     *              or better yet only accept constant V and J, two templated distributions,
     *              and a flag indicating if we should use the constants or the random distros
     */
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


