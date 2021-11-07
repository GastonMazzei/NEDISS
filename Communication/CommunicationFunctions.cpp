//
// Created by m4zz31 on 5/11/21.
//


#include "CommunicationFunctions.h"
#include "../Utils/global_standard_messages.h"

// Keep building MPI with https://www.boost.org/doc/libs/1_77_0/doc/html/mpi/tutorial.html

void irespond_value(ReferenceContainer &REF, double ix, int owner, std::list<MPI_Request>::iterator &R){
    auto vs = vertices(*REF.p_g);
    for (auto v = vs.first; v != vs.second; ++v){
        if (get(get(boost::vertex_index, *(REF.p_g)), *v) == ix){
            send_nonblocking(owner, *R, &(*REF.p_g)[*v].value, 3);
            return;
        }
    }
};


void send_nonblocking(int owner, MPI_Request &r, double *ix, int TAG){
    // (A) TAG = 1 means asking for a node's value by sending one int
    // (B) TAG = 2 means asking for an edge's value by sending two ints
    // (C) TAG = 3 means sending in return a node's value
    // (D) TAG = 4 means sending in return both an edge and node's value
    if ((TAG == 1) || (TAG == 3)) {
        MPI_Isend(ix,
                  1,//count
                  MPI_DOUBLE, // type
                  owner, // destination
                  TAG, MPI_COMM_WORLD, &r);
    } else if ((TAG == 2) || (TAG == 4)) {
        MPI_Isend(ix,
                  2,//count
                  MPI_DOUBLE, // type
                  owner, // destination
                  TAG, MPI_COMM_WORLD, &r);
        ++ix;
        MPI_Isend(ix,
                  2,//count
                  MPI_DOUBLE, // type
                  owner, // destination
                  TAG, MPI_COMM_WORLD, &r);
    }
};



void recv_nonblocking(int owner, MPI_Request &r, double *result, int TAG){
    // (A) TAG = 1 means getting a request for a node's value by recieving an int
    // (B) TAG = 2 means getting a request for an edge's value which require two ints
    // (C) TAG = 3 means getting in return a node's value
    // (D) TAG = 4 means getting in return both an edge and node's value
    if ((TAG == 1) || (TAG == 3)) {
        MPI_Irecv(result,
                  1,//count
                  MPI_DOUBLE, // type
                  owner, // destination
                  TAG, MPI_COMM_WORLD, &r);
    } else if ((TAG == 4) || (TAG == 2)) {
        MPI_Irecv(result,
                  2,//count
                  MPI_DOUBLE, // type
                  owner, // destination
                  TAG, MPI_COMM_WORLD, &r);
        ++result;
        MPI_Irecv(result,
                  2,//count
                  MPI_DOUBLE, // type
                  owner, // destination
                  TAG, MPI_COMM_WORLD, &r);
    }
};


void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g){
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

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g){
    //asks processor 'owner' for the 'v' and 'e' elements
};















