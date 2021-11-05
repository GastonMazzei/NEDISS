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
#include <set>
#include <boost/mpi/environment.hpp>


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

    // Instantiate the communication, parallel, integration and mapping helpers :-)
    //
    // !!! THEY SHOULD BE CONSUMED AS REFERENCES, NOT INSTANTIATED IN EACH CALL !!!
    //
    unsigned long NVtot = boost::num_vertices(g), NT;
    CommunicationHelper ComHelper(g);
    ParallelHelper ParHelper(ComHelper.NUM_THREADS, NVtot);
    IntegrationHelper IntHelper[NVtot];
    MappingHelper MapHelper(g);
    boost::mpi::environment env(boost::mpi::threading::funneled);



    std::set<long> CHECKED;
    NT = ComHelper.NUM_THREADS;
    auto vs = vertices(g);
    double temporalResult;

#pragma omp parallel firstprivate(NVtot, vs, NT)// Try me out! :-) Iterators are random access ones ;-)
{

    bool I_AM_NOT_MASTER = true; // this is mistrusting implementation of
                                // OpenMP standard where Master implies thread_id = 0 :-(
    long i=-1;
#pragma omp master
{
    mssleep(500);
    I_AM_NOT_MASTER = false;
};
    if (I_AM_NOT_MASTER) {

        long MY_THREAD = omp_get_thread_num();
        long N_THREADS = omp_get_num_threads();
        long MY_OFFSET, MY_LENGTH;
        long NLocals, M, rank;


        MY_OFFSET = (NVtot / (N_THREADS-1)) * (MY_THREAD-1) ; // (N_THREADS-1) instead of (N_THREADS)
        MY_LENGTH = (NVtot / (N_THREADS-1)); // because master will be entirely devoted to MPI ^.^
        if (MY_THREAD + 1 == N_THREADS) {
            MY_LENGTH += NLocals % N_THREADS;
        }

        i += MY_OFFSET;
        for (auto v = vs.first + MY_OFFSET; v != vs.first + (MY_OFFSET + MY_LENGTH - 1); ++v) {
            ++i;

            // Capture the central node's values
            IntHelper[i].centralValue = g[*v].value;
            IntHelper[i].centralParams = g[*v].params;
            // Get the neighbor and edge's values
            auto neighbors = boost::adjacent_vertices(*v, g);
            rank = in_degree(*v, g) + out_degree(*v, g);
            NLocals = neighbors.second - neighbors.first;
            M = rank - NLocals;


#pragma omp parallel firstprivate(i, M, neighbors, NLocals)
{
            long MY_THREAD_n = omp_get_thread_num();
            long N_THREADS_n = omp_get_num_threads();
            long MY_OFFSET_n, MY_LENGTH_n;
            MY_OFFSET_n = (NLocals / N_THREADS_n) * MY_THREAD_n;
            MY_LENGTH_n = (NLocals / N_THREADS_n);
            if (MY_THREAD_n + 1 == N_THREADS_n) {
                MY_LENGTH_n += NLocals % N_THREADS_n;
            }

            // debugging station ;-)
//            std::cout << "i: " << i << std::endl;
//            std::cout << "M: " << M << std::endl;
//            std::cout << "rank: " << rank << std::endl;
//            std::cout << "NLocals: " << NLocals << std::endl;
//            std::cout << "in_degree: " << in_degree(*v, g) << std::endl;
//            std::cout << "out_degree: " << out_degree(*v, g) << std::endl;
//            std::cout << "MY_THREAD: " << MY_THREAD << std::endl;
//            std::cout << "MY_OFFSET: " << MY_OFFSET << std::endl;
//            std::cout << "MY_LENGTH: " << MY_LENGTH << std::endl;
//            std::cout << "N_THREADS: " << N_THREADS << std::endl;
//            std::cout << "MY_THREAD_n: " << MY_THREAD_n << std::endl;
//            std::cout << "MY_OFFSET_n: " << MY_OFFSET_n << std::endl;
//            std::cout << "MY_LENGTH_n: " << MY_LENGTH_n << std::endl;
//            std::cout << "N_THREADS_n: " << N_THREADS_n << std::endl;


            if (M == 0) { // THE LEAST PROBABLE CASE! :-)
                // Skip directly to a Protocol to populate locally
                for (auto n = neighbors.first + MY_OFFSET_n;
                    n != neighbors.first + MY_OFFSET_n + MY_LENGTH_n; ++n) {
                    ED local_e = edge(*v, *n, g).first;
                    VD local_v = *n;
                    ParHelper.data[i].ProcessLocally[MY_THREAD_n].push_back(
                            std::make_pair(std::as_const(local_v), std::as_const(local_e)));
                }
            } else { // THE MOST PROBABLE CASE! :-)

                for (auto n = neighbors.first + MY_OFFSET_n;
                            n != neighbors.first + MY_OFFSET_n + MY_LENGTH_n; ++n) {
                    auto local_e = edge(*v, *n, g).first;
                    auto local_v = *n;
                    if (get(MapHelper.Local, *n) == 1) { // Case (1): locally available elements ;-)
                        if (edge(*v, *n, g).second == 1) {
                            ParHelper.data[i].ProcessLocally[MY_THREAD_n].push_back(
                                    std::make_pair(std::as_const(local_v), std::as_const(local_e)));
                        } else { error_report("Push back mechanism for local nodes has failed"); };
                    } else { // Case (2.A): We 'see' this neighbor because it is connected
                        // via an edge that we own. Get the edge and record who is the other node's owner
                        ParHelper.data[i].Missing[MY_THREAD_n].push_back(
                                std::make_tuple(std::as_const(local_v),
                                                std::as_const(local_e),
                                                get(MapHelper.NodeOwner, *n)));
                    }
                }
                // Finally please contemplate the 3rd case:
                // we have in-edges that are not available
                // locally so we have neighbors that cant
                // be indexed by 'neighbors'
                auto in_edges = boost::in_edges(*v, g);
                for (auto e = in_edges.first; e != in_edges.second; ++e) {
                    auto local_e = *e;
                    auto local_v = boost::source(*e, g);
                    ParHelper.data[i].Missing[MY_THREAD_n].push_back(
                            std::make_tuple(std::as_const(local_v), std::as_const(local_e),
                                            get(MapHelper.EdgeOwner, *e)));
                }
            }
} // end of the nested parallelism! :-)
            std::cout << "End of the neighbor iteration! starting new lap :-)" << std::endl;
#pragma atomic
            CHECKED.insert(i); // Adding the index to the list of checked indexes ;-)

        } // end of the for :-)
    } // end of the NON_MASTER_ONLY section :-)
} // end of the parallel construct

        // ******************HOW THE INTEGRATION CONTINUES  ;-)*****************
        //
        // Perform the evolution and store the result in the central
        // node's temporal register
//        g[*v].temporal_register = solver.evolve(centralValue,
//                                             centralParams,
//                                             neighborValues,
//                                             edgeValues);
//        clear_vectors(centralParams, neighborValues, edgeValues);





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










    // We could block operations until all nodes have seen their neighbors ;-)
    // :-) there is no Queue so we must be strict.
    // Parallel BGL uses the BSP model, which we could enforce by replacing a barrier
    // with a synchronization directive.
    synchronize(g.process_group());


    // For all nodes,
    // central_value :=  temporal_register
    register_to_value(g);

    // Block operations until all nodes have been updated ;-)
    adsync_barrier<0>();
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H