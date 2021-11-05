//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_GRAPHFUNCTIONS_H
#define CPPPROJCT_GRAPHFUNCTIONS_H

#include "../Solvers/GeneralSolver.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/memory_management.h"
#include "../Utils/msleep.h"
#include "../Utils/error.h"
#include "GeneralGraph.h"
#include <mpi.h>



// processors will iterate their "N=num_vertex(g)" nodes,
// and in each iteration they will either
//          (A)  (volume term)
//               iterate through their "NLocal=rank" locally-available neighbors
//               uploading the NLocal results into a container that can be used
//               to compute the central node's value using the solver evolver :-)
//          (B)  (surface term)
//               iterate through their "Nlocal<rank" locally-available neighbors,
//               while collecting "Nnonlocal" (what boost::graph::distributed call)
//               'vertex descriptors', 'edge descriptors', and the vertex's owner
//               in order to get involved in the following step which is requiring
//               a copy of this DynamicNode from the node's owner. Also, a new
//              iteration is performed over inward edges in order to collect the same
//              from them.
//
// It is in this context that we require keeping track of the surface (B) and volume (A) term's
// information using "N" (nodes in the process) containers that store a parallel cell object.

void register_to_value(Graph &g);

typedef Graph::vertex_descriptor VD;
typedef Graph::edge_descriptor ED;
typedef std::pair<VD, ED> InfoVecElem;
typedef std::tuple<VD, ED, int> PartialInfoVecElem;

struct ParallelCell{
    // Container of variable size that can keep track of different threads storing objects
    std::vector<std::list<InfoVecElem>> ProcessLocally;
    std::vector<std::list<PartialInfoVecElem>> Missing;
    explicit ParallelCell(){};
};

struct ParallelHelper{
    std::vector<ParallelCell> data; // initialized to NNodes x NT
    explicit ParallelHelper(int NT, unsigned long NNodes);
};

struct IntegrationHelper{
    double centralValue;
    std::vector<double> centralParams;
    std::vector<double> edgeValues;
    std::vector<double> neighborValues;
};

struct CommunicationHelper{
    // _NUM_THREADS can be captured from
    std::vector<int> WORLD_RANK, WORLD_SIZE, MY_NUM;
    int NUM_THREADS;
    boost::mpi::environment ENV;
    boost::mpi::communicator WORLD;
    explicit CommunicationHelper(Graph &g);
};

struct MappingHelper{
    // Relevant maps that should be passed as an instance
    // of a struct ;-)
    EdgeOwnerMap EdgeOwner;
    OwnerMap NodeOwner;
    LocalVertexMap Local;
    GlobalVertexMap Global;
    explicit MappingHelper(Graph &g): EdgeOwner(get(boost::edge_owner, g)),
                                     NodeOwner(get(boost::vertex_owner, g)),
                                     Local(get(boost::vertex_local, g)),
                                     Global(get(boost::vertex_global, g)){};

};




// Fully declared template
template<typename DIFFEQ, typename SOLVER> // e.g. DIFFEQ = NoiselessKuramoto, SOLVER = EulerSolver<NoiselessKuramoto>
void single_evolution(Graph &g, GeneralSolver<DIFFEQ,SOLVER> &solver){
    // There is an evolution operator that should compute the next
    // value and store it in the object's temporal register

    // said evolution operator should take (vector1, vector2, ...)
    //              where:
    //                  vector1:    <the value of the edges>
    //                  vector2:    <the value of the node at the other end of the edges>
    //                  vector3:    central node's parameters
    //                  value:      central node's value
    //                  ??MATRIX??: the parameter of each of the nodes, which should not be necessary
    //                              as it could be encoded in the edge.

    // TEST1: check that edges and nodes are correctly indexed, i.e.
    // (1)-1-(x)-2-(2) gets a vector1 and vector2 that is <1,2> and <1,2>

    // Instantiate the communication, parallel and mapping helpers :-)
    //
    // !!! THEY SHOULD BE CONSUMED AS REFERENCES, NOT INSTANTIATED IN EACH CALL !!!
    //
    CommunicationHelper ComHelper(g);
    ParallelHelper ParHelper(ComHelper.NUM_THREADS, num_vertices(g));
    MappingHelper MapHelper(g);


    // Auxiliaries for computations
    unsigned int ctrl_ownr, nghb_ownr;
    long i=-1;
    long NLocals, M, rank;
    VD local_v[ComHelper.NUM_THREADS]; // dynamic vertex and edge descriptors
    ED local_e[ComHelper.NUM_THREADS];  // (respectively)
    auto vs = vertices(g);
    double temporalResult;

    //---------------USE-ME-IF-YOU-NEED-TO-DEBUG--;-)-----------------
    //std::cout << "So far so good, flag N " << DEBUG_FLAG << std::endl;
    //++DEBUG_FLAG;
    //----------------------------------------------------------------
    // Start the loop! :-)
    for (auto v = vs.first; v != vs.second; ++v) {
        i = v - vs.first;


        // Instantiate the Integration Helper from this node's perspective ;-)
        const unsigned long Nmax = boost::num_vertices(g);
        IntegrationHelper IntHelper;

        // Capture the central node's values
        IntHelper.centralValue = g[*v].value;
        IntHelper.centralParams = g[*v].params;


        // Get the neighbor and edge's values
        auto neighbors = boost::adjacent_vertices(*v, g);
        rank = in_degree(*v, g) + out_degree(*v, g);
        NLocals = neighbors.second - neighbors.first;
        M = rank - NLocals;
        // THERE ARE 2 POSSIBLE SCENARIOS:
        // (1) all the data is available locally, which is Nlocals == rank
        // (2) there are M missing neighbors, where M = rank - NLocals,

        // Addressing (2):
        // there are 2 types of missing values.
        // 2.A) Nodes that are absent, but the edge object is available.
        //      This requires using that we know who is the owner and asking
        //      for the respective node.
        //      u=source(*ep.first,g);
        //      v=target(*ep.first,g);
        //      edge(vertex_descriptor u, vertex_descriptor v,
        //              const adjacency_list& g);
        //
        // 2.B) Nodes that are absent, and we only know about the incident edge.
        //      We must get the edge's owner and ask him for both the node and
        //      the edge's information.
        //      THIS DEFINES THE NEED FOR AN ASK_FOR_BOTH CALL


        if (M == 0) { // THE LEAST PROBABLE CASE! :-)
            // Skip directly to a Protocol to populate locally
            std::cout << "Implementing 'ALL-LOCAL' shortcut :-)" << std::endl;
#pragma omp parallel // Try me out! :-) Iterators are random access ones ;-)
{
                // ***[beginning][parallel]***
                // balancing workload among threads :O
                long MY_THREAD = omp_get_thread_num();
                long N_THREADS = omp_get_num_threads();
                long MY_OFFSET, MY_LENGTH;
                MY_OFFSET = (NLocals / N_THREADS) * MY_THREAD;
                MY_LENGTH = (NLocals / N_THREADS);
                if (MY_THREAD + 1 == N_THREADS) {
                    MY_LENGTH += NLocals % N_THREADS;
                }
                std::cout << "Hey there! I'm a thread of a parallel region A telling you that I'm thread N " << MY_THREAD << std::endl;
                for (auto n = neighbors.first + MY_OFFSET; n != neighbors.first + MY_OFFSET + MY_LENGTH; ++n) { // parallel
                    local_e[MY_THREAD] = edge(*v, *n, g).first;
                    local_v[MY_THREAD] = *n;//
                    std::cout << "About to write... " << std::endl;
                    ParHelper.data[i].ProcessLocally[MY_THREAD].push_back(std::make_pair(std::as_const(local_v[MY_THREAD]), std::as_const(local_e[MY_THREAD])));
                    std::cout << "...SUCCESS! " << std::endl;
                }
}
            std::cout << "Ended the 'ALL-LOCAL' shortcut :-)" << std::endl;
        } else { // THE MOST PROBABLE CASE! :-)
#pragma omp parallel // Try me out! :-) Iterators are random access ones ;-)
{
                // ***[beginning][parallel]***
                long MY_THREAD = omp_get_thread_num();
                long N_THREADS = omp_get_num_threads();
                long MY_OFFSET, MY_LENGTH;
                MY_OFFSET = (NLocals / N_THREADS) * MY_THREAD;
                MY_LENGTH = (NLocals / N_THREADS);
                if (MY_THREAD + 1 == N_THREADS) {
                    MY_LENGTH += NLocals % N_THREADS;
                }
                //ProcessLocally
                std::cout << "Hey there! I'm a thread of a parallel region B telling you that I'm thread N " << MY_THREAD << std::endl;
                // ***[end][parallel]***
                for (auto n = neighbors.first + MY_OFFSET; n != neighbors.first + MY_OFFSET + MY_LENGTH; ++n) {
                    if (get(MapHelper.Local, *n) == 1) { // Case (1): locally available elements ;-)
                        if (edge(*v, *n, g).second == 1) {
                            local_e[MY_THREAD] = edge(*v, *n, g).first;
                            local_v[MY_THREAD] = *n;//       "MY_OFFSET" here is the only parallel trait ;-) remove for sequential op
                        } else { error_report("Push back mechanism for local nodes has failed"); };
                    } else { // Case (2.A): We 'see' this neighbor because it is connected
                        // via an edge that we own. Get the edge and record who is the other node's owner
                        local_e[MY_THREAD] = edge(*v, *n, g).first;
                        local_v[MY_THREAD] = *n;
                        ParHelper.data[i].Missing[MY_THREAD].push_back(std::make_tuple(std::as_const(local_v[MY_THREAD]), std::as_const(local_e[MY_THREAD]),get(MapHelper.NodeOwner, *n)));
                    }
                }
            // Finally please contemplate the 3rd case:
            // we have in-edges that are not available
            // locally so we have neighbors that cant
            // be indexed by 'neighbors'
            auto in_edges = boost::in_edges(*v, g);
            for (auto e = in_edges.first; e != in_edges.second; ++e) {
                local_e[MY_THREAD] = *e;
                local_v[MY_THREAD] = boost::source(*e, g);
                ParHelper.data[i].Missing[MY_THREAD].push_back(std::make_tuple(std::as_const(local_v[MY_THREAD]), std::as_const(local_e[MY_THREAD]),get(MapHelper.EdgeOwner, *e)));
            }
} // end of parallel section ;-)
        }

        adsync_barrier<70>();
        std::cout << "End of the neighbor iteration! starting new lap :-)" << std::endl;
        adsync_barrier<30>();
    }


// **************************SOME COMMUNICATION WILL BE REQUIRED ;-)
//            if (world.rank() == 0) {
//                world.send(1, 0, std::string("Hello"));
//                std::string msg;
//                world.recv(1, 1, msg);
//                std::cout << msg << "!" << std::endl;
//            } else {
//                std::string msg;
//                world.recv(0, 0, msg);
//                std::cout << msg << ", ";
//                std::cout.flush();
//                world.send(0, 1, std::string("world"));
//            }

//            if ((world_rank == get(NodeOwner, *n)) && (world_rank != get(NodeOwner, *v))) {
////                MPI_Send(&g[*n].value, 1, MPI_DOUBLE, get(NodeOwner, *v), 0, MPI_COMM_WORLD);
//                world.send(get(NodeOwner, *v), 0, double(g[*n].value));
//            } else if ((world_rank == get(NodeOwner, *v)) && (world_rank != get(NodeOwner, *n))) {
//                std::cout << "I am process " << MY_NUM << " and I am trying to access a node I dont own "<< std::endl;
////                MPI_Recv(&temporalResult, 1, MPI_DOUBLE, get(NodeOwner, *n), 0, MPI_COMM_WORLD,
////                         MPI_STATUS_IGNORE);
//                world.recv(get(NodeOwner, *n), 0, temporalResult);
//                std::cout << "I recieved " << temporalResult << std::endl;
//                std::cout.flush();
//            }
    //adsync_barrier<20>();
    //}





        // ******************HOW THE INTEGRATION CONTINUES  ;-)*****************
        //
        // Perform the evolution and store the result in the central
        // node's temporal register
//        g[*v].temporal_register = solver.evolve(centralValue,
//                                             centralParams,
//                                             neighborValues,
//                                             edgeValues);

        // Empty the vectors using a function from utils :-)
        //clear_vectors(centralParams, neighborValues, edgeValues);

    //}








    // We could block operations until all nodes have seen their neighbors ;-)
    // :-) there is no Queue so we must be strict.
    // Parallel BGL uses the BSP model, which we could enforce by replacing a barrier
    // with a synchronization directive.
    //synchronize(g.process_group());


    // For all nodes,
    // central_value :=  temporal_register
    register_to_value(g);

    // Block operations until all nodes have been updated ;-)
    adsync_barrier<0>();
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H