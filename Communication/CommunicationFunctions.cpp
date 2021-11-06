//
// Created by m4zz31 on 5/11/21.
//


#include "CommunicationFunctions.h"
#include "../Utils/global_standard_messages.h"

// Keep building MPI with https://www.boost.org/doc/libs/1_77_0/doc/html/mpi/tutorial.html

void ask_for_node(int owner, VD &v, CommunicationHelper &H, int ix, Graph &g){
    // asks processor 'owner' for the 'v' element
    // while trying to dispatch as much requests
    // as possible ;-)
    //
    // This should be templated in order to use hyper-params comfortably,
    // e.g. batched-ops
    //
    // This could include a total pend in queue by reference and participate in
    // GetAllMsgs and GetOneMsg concurently...

    // Call "Nprocs - 1" probes with TAG=1 as in looking for questions
    int unread = H.WORLD_SIZE[0]-1;
    std::vector<int> flag(H.WORLD_SIZE[0],0);
    for (int i=0; i<H.WORLD_SIZE[0]; i++) {
        if  (i != H.WORLD_RANK[0]) {
            MPI_Iprobe(i, 1, MPI_COMM_WORLD, &flag[i], MPI_STATUS_IGNORE);
        }
    }

    // Send my specific request with TAG = 1 as it is a question
    MPI_Request r;
    int r_stat_1 = 0, r_stat_2 = 0;
    MPI_Isend(&ix,
              1,//count
              MPI_DOUBLE, // type
              owner, // destination
              1, MPI_COMM_WORLD, &r);

    // Loop over the probes waiting for some potential question with TAG=1
    while ((unread!=0) && (!r_stat_2)){
        // one lap checking for requests and responding them
        for (int i=0; i<H.WORLD_SIZE[0]; i++) {
            if  ((i != H.WORLD_RANK[0]) && flag[i]){
                double value;
                MPI_Recv(&value, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //---
                // respond accordingly!
                //--- USE TAG 0, FOR ANSWERS ONLY
                unread--;
            }
            if (!r_stat_1) {
                MPI_Request_get_status(r, &r_stat_1, MPI_STATUS_IGNORE);
                if (r_stat_1) {
                    double result;
                    r = MPI_Request();
                    //     Now we use TAG 0 which is the tag for answers :-)
                    MPI_Irecv(&result, 1, MPI_DOUBLE, owner, 0, MPI_COMM_WORLD, &r);
                }
            }
            if ((r_stat_1) && (!r_stat_2)){MPI_Request_get_status(r, &r_stat_2, MPI_STATUS_IGNORE);}
        }
    }
};

void ask_for_node_and_vertex(int owner, VD &v, ED &e, CommunicationHelper &H, int ix, Graph &g){
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
            ask_for_node(owner, new_v, H, their_vix, g);
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
            ask_for_node_and_vertex(owner, new_v, new_e, H, their_vix, g);
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