//
// Created by m4zz31 on 3/11/21.
//

#include "GraphFunctions.h"


void register_to_value(Graph &g){
    auto vs = vertices(g);
    for (auto v = vs.first; v != vs.second; ++v) {
        g[*v].value = g[*v].temporal_register;
    }
};

CommunicationHelper::CommunicationHelper(Graph &g) {
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
        data[i].Missing.resize(NT);
    }
}

