//
// Created by m4zz31 on 5/11/21.
//


#include "CommunicationFunctions.h"
#include "../Utils/global_standard_messages.h"

// Keep building MPI with https://www.boost.org/doc/libs/1_77_0/doc/html/mpi/tutorial.html

void ask_for_node(int owner, VD &v, CommunicationHelper &H, int ix){
    //asks processor 'owner' for the 'v' element
    if (H.WORLD_RANK[0] == 0) {
        H.WORLD.send(1, 0, std::string("Hello"));
        std::string msg;
        H.WORLD.recv(1, 1, msg);
        std::cout << msg << "!" << std::endl;
    } else {
        std::string msg;
        H.WORLD.recv(0, 0, msg);
        std::cout << msg << ", ";
        std::cout.flush();
        H.WORLD.send(0, 1, std::string("world"));
    }
};

void ask_for_node_and_vertex(int owner, VD &v, ED &e, CommunicationHelper &H, int ix){
    //asks processor 'owner' for the 'v' and 'e' elements


};

void GetOneMsg(int ix,
               CommunicationHelper &H,
               Graph &g,
               ParallelHelper &P,
               IntegrationHelper &I,
               std::queue<long> &C){
    std::list<InfoVecElem> Results;
    // Collect missing nodes
    for (auto & thread : P.data[ix].MissingA) {
        for (auto it =  thread.begin();
             it != thread.end(); it++){
            // retrieve data
            VD new_v;
            VD v = std::get<0>(*it);
            ED e = std::get<1>(*it);
            int our_vix = std::get<2>(*it);
            int their_vix = std::get<3>(*it);
            int owner = std::get<4>(*it);
            // ask for it, recieve it, and store it.
            ask_for_node(owner, new_v, H, their_vix);
#pragma atomic
            I[ix].ResultsPendProcess.emplace_back(new_v, std::as_const(e),our_vix, their_vix);
        }
    }
    // Collect missing nodes+edges
    for (auto & thread : P.data[ix].MissingB) {
        for (auto it =  thread.begin();
             it != thread.end(); it++){
            // retrieve data
            VD new_v;
            ED new_e;
            VD v = std::get<0>(*it);
            ED e = std::get<1>(*it);
            int our_vix = std::get<2>(*it);
            int their_vix = std::get<3>(*it);
            int owner = std::get<4>(*it);
            // ask for it, recieve it, and store it.
            ask_for_node_and_vertex(owner, new_v, new_e, H, their_vix);
#pragma atomic
            I[ix].ResultsPendProcess.emplace_back(new_v, std::as_const(e),our_vix, their_vix);
        }
    }
}



void GetAllMsgs(int NNodes,
               CommunicationHelper &H,
               Graph &g,
               ParallelHelper &P,
               IntegrationHelper &I,
               std::queue<long> &C){

    int READ = 1;
    long ix;
    const int TIME_COUNTER_TOLERANCE = 100; // 20 times delta_time
    int time_counter = 1, delta_time = 25;

    while (READ < NNodes) {
        while (C.empty()){
            if (time_counter > TIME_COUNTER_TOLERANCE){
                error_report(queue_timeout);
            }
            mssleep(delta_time);
            time_counter += 1;
        }
#pragma atomic
        ix = C.front();
        std::cout << "Queue was not empty! we got index: " << ix << std::endl;
#pragma atomic
        C.pop();
        std::cout << "About to iterate" << std::endl;
        GetOneMsg(ix,H,g,P,I,C);
        std::cout << "Ended!" << std::endl;
        time_counter = 1;
        ++READ;
    }
};
//
//while (condition)
//{
//MPI_Iprobe(...,&flag,...);
//if (flag)
//{
//MPI_Recv(...);
//...
//}
//// Do something, e.g. background tasks
//}