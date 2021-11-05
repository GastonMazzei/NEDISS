//
// Created by m4zz31 on 5/11/21.
//

#ifndef CPPPROJCT_HELPERCLASSES_H
#define CPPPROJCT_HELPERCLASSES_H

#include "../GraphClasses/GeneralGraph.h"


typedef Graph::vertex_descriptor VD;
typedef Graph::edge_descriptor ED;
typedef std::tuple<VD, ED, int, int> InfoVecElem;
typedef std::tuple<VD, ED, int, int, int> PartialInfoVecElem;

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

struct IntegrationCell{
    double centralValue;
    std::vector<double> centralParams;
    std::list<InfoVecElem> ResultsPendProcess;
    std::vector<double> edgeValues;
    std::vector<double> neighborValues;
};

typedef std::vector<IntegrationCell> IntegrationHelper;



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

struct OpenMPHelper{
    long MY_THREAD_n, N_THREADS_n, MY_OFFSET_n, MY_LENGTH_n;
    OpenMPHelper(long NLocals, int i);
    OpenMPHelper(long NLocals, int i, long N_THREADS, long MY_THREAD);
};



#endif //CPPPROJCT_HELPERCLASSES_H
