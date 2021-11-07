//
// Created by m4zz31 on 5/11/21.
//
#include "HelperClasses.h"

CommunicationHelper::CommunicationHelper(Graph &g) {
    // set the threading lvl: https://www.boost.org/doc/libs/1_68_0/doc/html/boost/mpi/threading/level.html
    // (our version is 1.77)
    //    boost::mpi::environment env(boost::mpi::threading::funneled);
    boost::mpi::environment env(boost::mpi::threading::serialized);

    int MY_NUM_v, WORLD_RANK_v, WORLD_SIZE_v;
    NUM_THREADS = std::stoi(std::getenv("OMP_THREAD_LIMIT"));
    MY_NUM_v = process_id(g.process_group());
    MPI_Comm_rank(MPI_COMM_WORLD, &WORLD_RANK_v);
    MPI_Comm_size(MPI_COMM_WORLD, &WORLD_SIZE_v);
    for (int i = 0; i < NUM_THREADS; i++) {
        WORLD_RANK.push_back(WORLD_RANK_v);
        WORLD_SIZE.push_back(WORLD_SIZE_v);
        MY_NUM.push_back(MY_NUM_v);
    }
}

ParallelHelper::ParallelHelper(int NT, unsigned long NNodes){
    data.resize(NNodes);
    for (int i=0; i<NNodes; i++){
        data[i].ProcessLocally.resize(NT);
        data[i].MissingA.resize(NT);
        data[i].MissingB.resize(NT);
    }
}


OpenMPHelper::OpenMPHelper(long NLocals, int i){
    MY_THREAD_n= omp_get_thread_num();
    N_THREADS_n = omp_get_num_threads();
    MY_OFFSET_n = (NLocals / (N_THREADS_n-i)) * (MY_THREAD_n-i);
    MY_LENGTH_n = (NLocals / (N_THREADS_n-i));
    if (MY_THREAD_n + 1 == N_THREADS_n) {
        MY_LENGTH_n += NLocals % (N_THREADS_n-i);
    }
}


OpenMPHelper::OpenMPHelper(long NLocals, int i, long N_THREADS, long MY_THREAD){
    MY_THREAD_n= MY_THREAD;
    N_THREADS_n = N_THREADS;
    MY_OFFSET_n = (NLocals / (N_THREADS_n-i)) * (MY_THREAD_n-i);
    MY_LENGTH_n = (NLocals / (N_THREADS_n-i));
    if (MY_THREAD_n + 1 == N_THREADS_n) {
        MY_LENGTH_n += NLocals % (N_THREADS_n-i);
    }
}


void IntegrationCell::build(Graph &g, VD v, MappingHelper &Map,
                            unsigned long &NOwned,
                            unsigned long &rank,
                            unsigned long &NLocals,
                            unsigned long &M){
    auto neighbors = boost::adjacent_vertices(v, g);
    auto in_edges = boost::in_edges(v, g);
    NOwned = 0;
    rank = in_degree(v, g) + out_degree(v, g);
    NLocals = neighbors.second - neighbors.first;
    M = rank - NLocals;
    for (auto n = neighbors.first; n != neighbors.second; ++n) {
        if (get(Map.Local, *n) == 1){
            ++NOwned;
        }
    }
    neighborValues.resize(rank);
    edgeValues.resize(rank);
}

ReferenceContainer::ReferenceContainer(ParallelHelper &ParHelper,
                                       CommunicationHelper &ComHelper,
                                       Graph & g,
                                       std::queue<long> & CHECKED,
                                       std::queue<long> & READY_FOR_INTEGRATION,
                                       IntegrationHelper & IntHelper,
                                       int & TOT,
                                       int & PENDING_INT){
    p_ParHelper = &ParHelper;
    p_ComHelper = &ComHelper;
    p_PENDING_INT = &PENDING_INT;
    p_g = &g;
    p_CHECKED = &CHECKED;
    p_READY_FOR_INTEGRATION = &READY_FOR_INTEGRATION;
    p_IntHelper = reinterpret_cast<IntegrationHelper *>(&IntHelper);
    p_TOT = &TOT;
}