//
// Created by m4zz31 on 5/11/21.
//

#ifndef CPPPROJCT_HELPERCLASSES_H
#define CPPPROJCT_HELPERCLASSES_H

#include "../GraphClasses/GeneralGraph.h"


typedef Graph::vertex_descriptor VD;
typedef Graph::edge_descriptor ED;//                        (UID is [Number of total nodes] * P_id  + local_index)
typedef std::tuple<double, double, unsigned long> InfoVecElem; // DynamicNode.value, DynamicEdge.value, UID
typedef std::tuple<double, int, int> PartialInfoVecElem; // DynamicNode.value Owner Index

struct ParallelCell{
    // Container of variable size that can keep track of different threads storing objects
    std::vector<std::list<InfoVecElem>> ProcessLocally;
    std::vector<std::list<PartialInfoVecElem>> MissingA;
    std::vector<std::list<PartialInfoVecElem>> MissingB;
    explicit ParallelCell(){};
};


struct ParallelHelper{
    std::vector<ParallelCell> data; // initialized to NNodes x NT
    explicit ParallelHelper(int NT, unsigned long NNodes);
};

struct LayeredSolverCell{
    std::vector<double> RK1;
    std::vector<double> RK2;
    std::vector<double> RK3;
    std::vector<double> RK4;
    std::vector<bool> RK1_status;
    std::vector<bool> RK2_status;
    std::vector<bool> RK3_status;
    std::vector<bool> RK4_status;
    LayeredSolverCell(int N);
    LayeredSolverCell(){};
};

struct LayeredSolverHelper : public  LayeredSolverCell {
public:
    bool built = false;
    std::vector<LayeredSolverCell> data;
    LayeredSolverHelper(int N): data(N){};
    void buildForRank(long ix, long rank);
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



struct IntegrationCell{
    double centralValue;
    std::vector<double> centralParams;
    std::list<InfoVecElem> ResultsPendProcess;
    std::vector<double> edgeValues;
    std::vector<double> neighborValues;
    void build(Graph &g, VD v, MappingHelper &Map,
               unsigned long &NOwned,
               unsigned long &rank,
               unsigned long &NLocals,
               unsigned long &M);
};


typedef std::vector<IntegrationCell> IntegrationHelper;


struct ReferenceContainer {
    double placeholder = 0.;
    ParallelHelper * p_ParHelper;
    CommunicationHelper * p_ComHelper;
    IntegrationHelper * p_IntHelper;
    MappingHelper * p_MapHelper;
    LayeredSolverHelper * p_LayHelper;
    Graph * p_g;
    int * p_TOT;
    bool * p_keepResponding;
    int * p_PENDING_INT;
    std::pair<std::queue<long>, std::queue<unsigned long>> * p_CHECKED;
    std::pair<std::queue<long>, std::queue<unsigned long>> * p_READY_FOR_INTEGRATION;
    ReferenceContainer(ParallelHelper &ParHelper,
                       CommunicationHelper &ComHelper,
                       Graph & g,
                       std::pair<std::queue<long>, std::queue<unsigned long>> & CHECKED,
                       std::pair<std::queue<long>, std::queue<unsigned long>> & READY_FOR_INTEGRATION,
                       IntegrationHelper & IntHelper,
                       int & TOT,
                       int & PENDING_INT,
                       MappingHelper & MapHelper,
                       LayeredSolverHelper & LayHelper,
                       bool keepResponding);
};


struct OpenMPHelper{
    long MY_THREAD_n, N_THREADS_n, MY_OFFSET_n, MY_LENGTH_n;
    OpenMPHelper(long NLocals, int i);
    OpenMPHelper(long NLocals, int i, long N_THREADS, long MY_THREAD);
};



#endif //CPPPROJCT_HELPERCLASSES_H
