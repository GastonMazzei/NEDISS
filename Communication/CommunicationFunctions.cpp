//
// Created by m4zz31 on 5/11/21.
//


#include "CommunicationFunctions.h"
#include "../Utils/global_standard_messages.h"

void build_answer_edges(double * answer, ReferenceContainer &REF, double * ix, int owner, int MyNProc){
    // Iterate through nodes
    bool found = false;
    auto vs = vertices(*REF.p_g);
    auto v = vs.first;
    while ((v != vs.second) && (!found)){
        if (get(REF.p_MapHelper->NodeOwner,*v)==MyNProc) {

            // If I am the owner, test if this node has the index I'm looking for
            if (get(get(boost::vertex_index, *(REF.p_g)), *v) == (int) *(ix+1)) {

                // If this is the node this code was aimed to, iterate through neighbrors
                auto neighbors = boost::adjacent_vertices(*v, *REF.p_g);
                auto n = neighbors.first;
                while ((n != neighbors.second) && (!found)) {

                    //Iterate through neighbors looking for the correct edge ;-)
                    if (get(REF.p_MapHelper->NodeOwner, *n) == owner) {

                        // If this neighbor is owned by the requester
                        if (get(get(boost::vertex_index, *(REF.p_g)), *n) == (int) *ix) {

                            // If its index is the one the requester refered to
                            found = true;
                            PRINTF_DBG("Found the requested node attached to the requested edge!\n");
                            auto e = edge(*v, *n, *(REF.p_g));
                            if (e.second != 1) {
                                printf("[FATAL] Requested edge did not exist!!!");
                                exit(1);
                            }
                            int TAG = OFFSET + (int) *ix;
                            *answer = (*(REF.p_g))[*v].value;
                            *(answer + 1) = (*(REF.p_g))[e.first].value;
                        }
                    }
                    ++n;
                }
            }
        }
        ++v;
    }
    if (!found) {
        printf("CRITICAL: node+edge not found!\n");
        std::cout << std::flush;
        PRINTF_DBG("[CRITICAL] THE NODE WAS NOT FOUND\n");
        std::cout << std::flush;
        exit(1);
    }
}



void build_answer(double &answer, ReferenceContainer &REF, double ix, int owner, int MyNProc){
    auto vs = vertices(*REF.p_g);
    for (auto v = vs.first; v != vs.second; ++v){
        if (get(REF.p_MapHelper->NodeOwner,*v)==MyNProc) {
            if (get(get(boost::vertex_index, *(REF.p_g)), *v) == (int) ix){
                answer = (*REF.p_g)[*v].value;
                return;
            }
        }
    }
    PRINTF_DBG("[CRITICAL] THE NODE WAS NOT FOUND\n");std::cout<<std::flush;
    exit(1);
}




//void build_answer_edges(double * answer, ReferenceContainer &REF, double * ix, int owner, int MyNProc){
//
//    // Iterate through nodes
//    bool found = false;
//    auto vs = vertices(*REF.p_g);
//    auto v = vs.first;
//    while ((v != vs.second) && (!found)){
//        if (get(REF.p_MapHelper->NodeOwner,*v)==MyNProc) {
//
//            // If I am the owner, test if this node has the index I'm looking for
//            if (get(get(boost::vertex_index, *(REF.p_g)), *v) == *(ix+1)) {
//
//                // If this is the node this code was aimed to, iterate through neighbrors
//                auto neighbors = boost::adjacent_vertices(*v, *REF.p_g);
//                auto n = neighbors.first;
//                while ((n != neighbors.second) && (!found)) {
//
//                    //Iterate through neighbors looking for the correct edge ;-)
//                    if (get(REF.p_MapHelper->NodeOwner, *n) == owner) {
//
//                        // If this neighbor is owned by the requester
//                        if (get(get(boost::vertex_index, *(REF.p_g)), *v) == *ix) {
//
//                            // If its index is the one the requester refered to
//                            found = true;
//                            printf("Found the requested node attached to the requested edge!\n");
//                            auto e = edge(*v, *n, *(REF.p_g));
//                            if (e.second != 1) {
//                                printf("[FATAL] Requested edge did not exist!!!");
//                                exit(1);
//                            }
//                            int TAG = OFFSET + *ix;
//                            *answer = (*(REF.p_g))[*v].value;
//                            *(answer + 1) = (*(REF.p_g))[e.first].value;
//                            printf("Found it OK\n"); std::cout<<std::flush;
//                            return;
//                        }
//                    }
//                    ++n;
//                }
//            }
//        }
//        ++v;
//    }
//    printf("CRITICAL: not found!\n"); std::cout<<std::flush;
//    PRINTF_DBG("[CRITICAL] THE NODE WAS NOT FOUND\n");std::cout<<std::flush;
//    exit(1);
//}




void respond_value(ReferenceContainer &REF, double ix, int owner, int MyNProc){
    auto vs = vertices(*REF.p_g);
    for (auto v = vs.first; v != vs.second; ++v){
        if (get(REF.p_MapHelper->NodeOwner,*v)==MyNProc) {
            //printf("I am the owner! :-)");                        // NEW: added (int) in next line
            if (get(get(boost::vertex_index, *(REF.p_g)), *v) == (int) ix){
                printf("About to respond synchronically and blockingly\n"); std::cout<<std::flush;
                MPI_Ssend(&(*REF.p_g)[*v].value, 1, MPI_DOUBLE, owner, (int) ix, MPI_COMM_WORLD);
                PRINTF_DBG("ANSWERED ONE MESSAGE!\n");
                PRINTF_DBG("Answering %d with %d's value, which is %f\n",owner, (int) ix, (*REF.p_g)[*v].value);
                return;
            } else {
                PRINTF_DBG("I am not the owner!!\n");
                std::cout  << std::flush;
            }
        }
    }
    PRINTF_DBG("[CRITICAL] THE NODE WAS NOT FOUND\n");std::cout<<std::flush;
    exit(1);
};

void irespond_value(ReferenceContainer &REF, double ix, int owner, MPI_Request & R, int MyNProc){
    auto vs = vertices(*REF.p_g);
    for (auto v = vs.first; v != vs.second; ++v){
        if (get(REF.p_MapHelper->NodeOwner,*v)==MyNProc) {
            //printf("I am the owner! :-)");                        // NEW: added (int) in next line
            if (get(get(boost::vertex_index, *(REF.p_g)), *v) == (int) ix){
                //send_nonblocking(owner, R, REF.placeholder, (int) ix);
                send_nonblocking(owner, R, (*REF.p_g)[*v].value, (int) ix);
                PRINTF_DBG("ANSWERED ONE MESSAGE!\n");
                PRINTF_DBG("Answering %d with %d's value, which is %f\n",owner, (int) ix, (*REF.p_g)[*v].value);
                return;
        } else {
            PRINTF_DBG("I am not the owner!!\n");
            std::cout  << std::flush;
        }
        }
    }
    PRINTF_DBG("[CRITICAL] THE NODE WAS NOT FOUND\n");std::cout<<std::flush;
    exit(1);
};










void send_nonblocking2(int owner, MPI_Request &r, double &ix, int TAG){
    int answer;
    answer = MPI_Isend(&ix,
                       2,//count
                       MPI_DOUBLE, // type
                       owner, // destination
                       TAG, MPI_COMM_WORLD, &r);
    PRINTF_DBG("Just sent a message with status %d\n", answer);
    std::cout << std::flush;
    if (answer!=0){
        PRINTF_DBG("Answer was a bad status, recursive reset\n");
        PRINTF_DBG("bypassed into direct exit: error code was %d\n",answer);
        std::cout << std::flush;
        exit(answer);
        send_nonblocking2(owner, r, ix, TAG);
    }
};






void send_nonblocking(int owner, MPI_Request &r, double &ix, int TAG){
    int answer;
    answer = MPI_Isend(&ix,
              1,//count
              MPI_DOUBLE, // type
              owner, // destination
              TAG, MPI_COMM_WORLD, &r);
    PRINTF_DBG("Just sent a message with status %d\n", answer);
    std::cout << std::flush;
    if (answer!=0){
        PRINTF_DBG("Answer was a bad status, recursive reset\n");
        PRINTF_DBG("bypassed into direct exit: error code was %d\n",answer);
        std::cout << std::flush;
        exit(answer);
        send_nonblocking(owner, r, ix, TAG);
    }
};



void recv_blocking(int owner, MPI_Request &r, double &result, int TAG){
    int answer;
    answer = MPI_Recv(&result,
                       1,//count
                       MPI_DOUBLE, // type
                       owner, // destination
                       TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (answer!=0) {
        PRINTF_DBG("Recieve was a bad status, recursive reset\n");
        recv_blocking(owner, r, result, TAG);
    }
}

void send_blocking(int owner, MPI_Request &r, double &ix, int TAG) {
    int answer;
    answer = MPI_Send(&ix,
                       1,//count
                       MPI_DOUBLE, // type
                       owner, // destination
                       TAG, MPI_COMM_WORLD);
    PRINTF_DBG("Just sent a (blocking) message with status %d\n", answer);
    std::cout << std::flush;
    if (answer != 0) {
        PRINTF_DBG("Answer was a bad status, recursive reset\n");
        PRINTF_DBG("bypassed into direct exit: error code was %d\n", answer);
        std::cout << std::flush;
        exit(answer);
        send_nonblocking(owner, r, ix, TAG);
    }
}

void recv_nonblocking2(int owner, MPI_Request &r, double & result, int TAG){
    int answer;
    answer = MPI_Irecv(&result,
                       2,//count
                       MPI_DOUBLE, // type
                       owner, // destination
                       TAG, MPI_COMM_WORLD, &r);

    if (answer!=0){
        PRINTF_DBG("Recieve was a bad status, recursive reset\n");
        recv_nonblocking2(owner, r, result, TAG);
    }
};


void recv_nonblocking(int owner, MPI_Request &r, double &result, int TAG){
    int answer;
    answer = MPI_Irecv(&result,
                  1,//count
                  MPI_DOUBLE, // type
                  owner, // destination
                  TAG, MPI_COMM_WORLD, &r);

    if (answer!=0){
        PRINTF_DBG("Recieve was a bad status, recursive reset\n");
        recv_nonblocking(owner, r, result, TAG);
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
        if  (i != H.WORLD_RANK[0]) { // // TAAAAG IS INCORRECT TAG IS INCORRECT !!! WARNING
            MPI_Iprobe(i, 1, MPI_COMM_WORLD, &flag[i], MPI_STATUS_IGNORE);// TAG 1 IS INCORRECT
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
                PRINTF_DBG("[HANG ALERT] About to use a blocking recieve in ask_for_node");
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

void destroyRequest(MPI_Request &R, int &NERR){
    MPI_Status S;
    int flag=0;
    MPI_Cancel(&R);
    MPI_Wait(&R, &S);
    MPI_Test_cancelled(&S, &flag);
    if (!flag){
        PRINTF_DBG("failed to cancel a request!  :-(");
        ++NERR;
    }
}


void destroyRequestWithoutCounter(MPI_Request &R){
    MPI_Status S;
    int flag=0;
    MPI_Cancel(&R);
    MPI_Wait(&R, &S);
    MPI_Test_cancelled(&S, &flag);
    if (!flag){
        PRINTF_DBG("\n\n\n\nfailed to cancel a request! flag was %d :-(\n\n\n\n",flag);
    }
}

void freeRequestWithoutCounter(MPI_Request &R){
    int statusFree;
    statusFree = MPI_Request_free(&R);
    if (statusFree != 0){
        PRINTF_DBG("\n\n\n\nfailed to FREE a request! status was %d :-(\n\n\n\n",statusFree);
    }
}

int destroyRequestReturnInteger(MPI_Request &R){
    MPI_Status S;
    int flag=0;
    MPI_Cancel(&R);
    MPI_Wait(&R, &S);
    MPI_Test_cancelled(&S, &flag);
    return flag;
}




void sendReqForTest(int MYPROC, int i){
    // Send request for missing info
    double vix = (double) MYPROC;
    double vval;
    int status_flagstatus=1;
    int issend_status = 1;
    int MAXTRIES = 20; // HYPERPARAM
    int counter = 0;
    int s_fs=0, s_fr=0;
    MPI_Request sReq, rReq;
    int owner, flagstatus;
    if (MYPROC==0){
        owner = 2;
    } else {
        owner = MYPROC -1;
    }


    issend_status = MPI_Issend(&vix, 1, MPI_DOUBLE, owner, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &sReq);
    while (issend_status != 0){
        issend_status = MPI_Issend(&vix, 1, MPI_DOUBLE, owner, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &sReq);
    }

    status_flagstatus = MPI_Test(&sReq, &flagstatus, MPI_STATUS_IGNORE);
    while ((!(status_flagstatus==0)) || (!(flagstatus==1))){
        status_flagstatus = MPI_Test(&sReq, &flagstatus, MPI_STATUS_IGNORE);
        mssleep(50); // doing something here :-)
    }

    PRINTF_DBG("I have supposedly sent this message to %d with no status_flagstatus nor flagstatus \n",
               owner);

    recv_nonblocking(owner,
                     rReq, // flag = (int) index :-)
                     vval, (int) vix);
    //mssleep(200);
    status_flagstatus = MPI_Test(&rReq, &flagstatus, MPI_STATUS_IGNORE);

    while ((!(status_flagstatus==0)) || (!(flagstatus==1))){
        counter++;
        PRINTF_DBG("Stuck at some special state B;-/ recieving %d and flag %d\n",status_flagstatus, flagstatus);
        status_flagstatus = MPI_Test(&rReq, &flagstatus, MPI_STATUS_IGNORE);
        //mssleep(50);
        if(counter >= MAXTRIES) {
            status_flagstatus = 0;
            flagstatus = 1;
        }
    }

    if (counter>=MAXTRIES){
        MPI_Cancel(&rReq);
        //MPI_Cancel(&sReq);
        PRINTF_DBG("RECURSIVE SOLUTION... the damned proc is %d\n",MYPROC);std::cout<<std::flush;
        sendReqForTest(MYPROC, i);
    }
    else {
        PRINTF_DBG("Sent and recieved correctly ;-)  (i=%d)\n", i);
        std::cout << std::flush;
    }


}










