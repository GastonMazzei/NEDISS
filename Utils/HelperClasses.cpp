//
// Created by m4zz31 on 5/11/21.
//
#include "HelperClasses.h"

CommunicationHelper::CommunicationHelper(Graph &g) {
    boost::mpi::environment env(boost::mpi::threading::funneled);
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
        MY_LENGTH_n += NLocals % N_THREADS_n;
    }
}


OpenMPHelper::OpenMPHelper(long NLocals, int i, long N_THREADS, long MY_THREAD){
    MY_THREAD_n= MY_THREAD;
    N_THREADS_n = N_THREADS;
    MY_OFFSET_n = (NLocals / (N_THREADS_n-i)) * (MY_THREAD_n-i);
    MY_LENGTH_n = (NLocals / (N_THREADS_n-i));
    if (MY_THREAD_n + 1 == N_THREADS_n) {
        MY_LENGTH_n += NLocals % N_THREADS_n;
    }
}