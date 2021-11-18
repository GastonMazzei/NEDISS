//
// Created by m4zz31 on 5/11/21.
//

#ifndef CPPPROJCT_COMMUNICATIONFUNCTIONS_H
#define CPPPROJCT_COMMUNICATIONFUNCTIONS_H
#include "../macros/macros.h"
#include <boost/serialization/string.hpp>
#include "../GraphClasses/GeneralGraph.h"
#include "../GraphClasses/GraphFunctions.h"
#include "../Utils/HelperClasses.h"
#include "../Utils/error.h"
#include "../Utils/global_standard_messages.h"
#include "../Utils/msleep.h"
#include <set>
#include <iterator>
#include <random>

void sendReqForTest(int MYPROC, int i);
int destroyRequestReturnInteger(MPI_Request &R);
void destroyRequest(MPI_Request &R, int &NERR);
void destroyRequestWithoutCounter(MPI_Request &R);
void freeRequestWithoutCounter(MPI_Request &R);
void send_nonblocking(int owner, MPI_Request &r, double &ix, int TAG);
void recv_nonblocking(int owner, MPI_Request &r, double &result, int TAG);
void recv_nonblocking2(int owner, MPI_Request &r, double & result, int TAG);
void irespond_value(ReferenceContainer &REF, double ix, int owner, MPI_Request & R, int MyNProc);
void respond_value(ReferenceContainer &REF, double ix, int owner, int MyNProc);
void irespond_value_edges(ReferenceContainer &REF, double *ix, int owner, MPI_Request & R, int MyNProc);
void send_nonblocking2(int owner, MPI_Request &r, double &ix, int TAG);
void build_answer(double &answer, ReferenceContainer &REF, double ix, int owner, int MyNProc);

// OLD ANSWER MESSAGES
//
//template<int DT, int TIMETOL, int BATCH>
//void answer_messages(ReferenceContainer &REF,int MYTHR) {
//
//    int MYPROC = REF.p_ComHelper->WORLD_RANK[MYTHR];
//    int NPROCS = REF.p_ComHelper->WORLD_SIZE[MYTHR];
//    std::vector<MPI_Request> R_tot;
//    std::vector<MPI_Request> R_send;
//    std::vector<double> ix;
//    int NDISPATCHED = 0;
//
//    // lay the probes for all the p rocs
//    int flagprobes[NPROCS];
//    int probe_status = 1;
//    int statusreq_status = 1;
//    for (int i = 0; i < NPROCS; i++) {
//        flagprobes[i] = 0;
//        if (i != MYPROC) {
//            probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//            while (probe_status!=0){
//                probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//            }
//        }
//    }
//
//    // initialize auxiliary variables
//    std::uniform_int_distribution<int> gen(0, NPROCS - 1);   // This three lines could be
//    unsigned int SEED = std::stoi(std::getenv("SEED"));   // encapsulated inside REF
//    std::mt19937 rng(SEED);                              // :-) and we make it 'more efficient' lol
//    int ticks = 0;
//    std::set<int> answered;
//    answered.insert(MYPROC);
//    //double ix;
//    int i = gen(rng); // initial processor to which check if we can respond to
//    if (i == MYPROC){
//        ++i;
//        if (i>=NPROCS){
//            i=0;
//        }
//    }
//    int statusFree = 0;
//    int NTOT = 0;
//    int status = 0;
//    int status_localreq = 0;
//    int NERR = 0;
//
//    for (int tk = 0; tk < TIMETOL; ++tk) {
//        for (int j = 0; j < BATCH; ++j) {
//            // Re-probe it :-)
//            probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//            while (probe_status!=0){
//                probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//            }
//
//            if (flagprobes[i] == 1) { // a message appears to be available! recieve it!
//
//                // Print the first three elements in the asynchronous probe :-)
//                PRINTF_DBG("|%d %d %d|\n",flagprobes[0] , flagprobes[1] , flagprobes[2]);std::cout<<std::flush;
//
//                // Try to recieve it
//                R_tot.push_back(MPI_Request());
//                ix.push_back((double) 0);
//                recv_nonblocking(i, R_tot[R_tot.size() - 1], ix[R_tot.size() - 1], VERTEXVAL_REQUEST_FLAG);
//
//                status_localreq = 0;
//                int localtimecounter = 0;
//                while ((status_localreq != 1) && (localtimecounter < TIMETOL)){
//                    statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
//                    while (statusreq_status != 0) {
//                        statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
//                    }
//                    localtimecounter++;
//                    //mssleep(DT);
//                }
//
//
//                if (status_localreq == 1) { // If we were first to capture the message, proceed.
//                    R_send.push_back(MPI_Request());
//                    PRINTF_DBG("We effectively captured a vertex info request :-)\n");
//                    irespond_value(REF, ix[R_tot.size() - 1], i, R_send[R_send.size() - 1], MYPROC);
//                    ++NTOT;
//                    PRINTF_DBG("We effectively answered asynchronously a vertex info request :-)\n");
//                    NDISPATCHED++;
//                } else {
//                    // MPI_Test didnt destroy the request, so explicitly free it
//                    statusFree = MPI_Request_free(&R_tot[R_tot.size() - 1]);
//                    PRINTF_DBG("We were faced with a probe that indicated an incoming message but we couldnt capture it :o\n");
//                }
//
//                // Reset vars.
//                status_localreq = 0;
//                flagprobes[i] = 0;
//            }
//
//
//
//            // For next iteration
//            ++i;
//            if (i == MYPROC) ++i;
//            if (i >= NPROCS) {
//                if (MYPROC!=0) {
//                    i = 0;
//                } else {
//                    i = 1;
//                }
//            }
//        }
//
//        // END OF THE BATCH... now we wait for all the sent requests. (In total they can be up to "Batch")
//        if (R_send.size() > 0) {
//
//            // WHY NOT JUST WAITALL
//
//            PRINTF_DBG("ENDing answer_messages. Caught (R_tot)=%d requests, Answered (R_send)=%d.\n",
//                       R_tot.size(),R_send.size());
//
////            int status_of_getstatus = 1;
////            for (int i = 0; i < R_send.size(); ++i) {
////                status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
////                while (status_of_getstatus != 0) {
////                    PRINTF_DBG("Failing to test a request, its returning %d\n", status_of_getstatus);
////                    status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
////                }
////
////                if (status_localreq == 0) {
////                    PRINTF_DBG("We are currently waiting for a send to be completed\n");
////                    MPI_Wait(&R_send[i], MPI_STATUS_IGNORE);
////                    PRINTF_DBG("Successfully waited the request to be completed\n");
////                }
////            }
//        }
//
//        std::vector<MPI_Request> R_tot;
//        std::vector<MPI_Request> R_send;
//        std::vector<double> ix;
//        //mssleep(DT);
//        ++ticks;
//    }
//    //printf("---answer_request heartbeat---\n");std::cout<<std::flush;
//};



template<int DT, int TIMETOL, int BATCH>
void answer_messages(ReferenceContainer &REF,int MYTHR) {

    int MYPROC = REF.p_ComHelper->WORLD_RANK[MYTHR];
    int NPROCS = REF.p_ComHelper->WORLD_SIZE[MYTHR];
    int NDISPATCHED = 0;
    int flag;
    int probe_status = 1;
    double ix;
    int statusreq_status = 1;
    double buffer;
    double answer;
    int ticks = 0;
    int statusFree = 0;
    int NTOT = 0;
    int NERR = 0;

    MPI_Request R;
    int TRIES=0;
    int t=0;
    bool firstlap=true;
    while (TRIES < 3){
        MPI_Status status;
        if ((flag == 1) || firstlap) {
            R = MPI_Request();
            MPI_Irecv(&buffer, 1, MPI_DOUBLE,
                      MPI_ANY_SOURCE,
                      VERTEXVAL_REQUEST_FLAG,
                      MPI_COMM_WORLD,
                      &R);
        }
        if (firstlap) firstlap = false;
        MPI_Request_get_status(R, &flag, &status);
        while ((flag != 1) && (t<TIMETOL)){
            MPI_Request_get_status(R, &flag, &status);
            ++t;
            mssleep(DT);
        }
        if (t>=TIMETOL) ++TRIES;
        t=0;
        if (flag == 1){
            build_answer(answer, REF, buffer, status.MPI_SOURCE, MYPROC);
            printf("About to send one!\n");std::cout<<std::flush;
            MPI_Ssend(&answer, 1, MPI_DOUBLE, status.MPI_SOURCE, (int) buffer, MPI_COMM_WORLD);
            printf("Correctly sent one!\n");std::cout<<std::flush;
        }
    }
    MPI_Cancel(&R);
    MPI_Status status;
    MPI_Request_get_status(R, &flag, &status);
    MPI_Test_cancelled(&status, &flag);
    if (flag!=1) {
        //MPI_Wait(&R, MPI_STATUS_IGNORE)
        printf("A request failed to be cancelled, we are assuming we recieved it!\n");std::cout<<std::flush;
        build_answer(answer, REF, buffer, status.MPI_SOURCE, MYPROC);
        MPI_Ssend(&answer, 1, MPI_INT, status.MPI_SOURCE, (int) buffer, MPI_COMM_WORLD);
    }
};



template<int DT, int TIMETOL, int BATCH>
void answer_messages_edges(ReferenceContainer &REF,int MYTHR) {
    // (1) recieve two numbers: my and the other's node index :K
    // (2) send through channel OFFSET+OWNER_INDEX  [vertexval, edgeval]

    int MYPROC = REF.p_ComHelper->WORLD_RANK[MYTHR];
    int NPROCS = REF.p_ComHelper->WORLD_SIZE[MYTHR];
    std::vector<MPI_Request> R_tot;
    std::vector<MPI_Request> R_send;
    int NDISPATCHED = 0;

    // lay the probes for all the p rocs
    int flagprobes[NPROCS];
    int probe_status = 1;
    int statusreq_status = 1;
    for (int i = 0; i < NPROCS; i++) {
        flagprobes[i] = 0;
        if (i != MYPROC) {
            probe_status = MPI_Iprobe(i, EDGEVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
            while (probe_status!=0){
                probe_status = MPI_Iprobe(i, EDGEVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
            }
        }
    }

    // initialize auxiliary variables
    std::uniform_int_distribution<int> gen(0, NPROCS - 1);   // This three lines could be
    unsigned int SEED = std::stoi(std::getenv("SEED"));   // encapsulated inside REF
    std::mt19937 rng(SEED);                              // :-) and we make it 'more efficient' lol
    int ticks = 0;
    std::set<int> answered;
    answered.insert(MYPROC);
    double ix[2];
    int i = gen(rng); // initial processor to which check if we can respond to
    if (i == MYPROC){
        ++i;
        if (i>=NPROCS){
            i=0;
        }
    }
    int statusFree = 0;
    int NTOT = 0;
    int status = 0;
    int status_localreq = 0;
    int NERR = 0;

    for (int tk = 0; tk < TIMETOL; ++tk) {
        for (int j = 0; j < BATCH; ++j) {
            if (flagprobes[i] == 1) { // a message appears to be available! recieve it!

                // Print the first three elements in the asynchronous probe :-)
                PRINTF_DBG("|%d %d %d|\n",flagprobes[0] , flagprobes[1] , flagprobes[2]);std::cout<<std::flush;

                // Try to recieve it
                R_tot.push_back(MPI_Request());
                recv_nonblocking2(i, R_tot[R_tot.size() - 1], ix[0], EDGEVAL_REQUEST_FLAG);

                status_localreq = 0;
                int localtimecounter = 0;
                while ((status_localreq != 1) && (localtimecounter < TIMETOL)){
                    statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
                    while (statusreq_status != 0) {
                        statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
                    }
                    localtimecounter++;
                    //mssleep(DT);
                }


                if (status_localreq == 1) { // If we were first to capture the message, proceed.
                    R_send.push_back(MPI_Request());
                    PRINTF_DBG("We effectively captured a vertex info request :-)\n");
                    irespond_value_edges(REF, &ix[0], i, R_send[R_send.size() - 1], MYPROC);
                    ++NTOT;
                    PRINTF_DBG("We effectively answered asynchronously a vertex info request :-)\n");
                    NDISPATCHED++;
                } else {
                    // MPI_Test didnt destroy the request, so explicitly free it
                    statusFree = MPI_Request_free(&R_tot[R_tot.size() - 1]);
                    PRINTF_DBG("We were faced with a probe that indicated an incoming message but we couldnt capture it :o\n");
                }

                // Reset vars.
                status_localreq = 0;
                flagprobes[i] = 0;
            }

            // Re-probe it :-)
            probe_status = MPI_Iprobe(i, EDGEVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
            while (probe_status!=0){
                probe_status = MPI_Iprobe(i, EDGEVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
            }

            // For next iteration
            ++i;
            if (i == MYPROC) ++i;
            if (i >= NPROCS) {
                if (MYPROC!=0) {
                    i = 0;
                } else {
                    i = 1;
                }
            }
        }

        // END OF THE BATCH... now we wait for all the sent requests. (In total they can be up to "Batch")
        if (R_send.size() > 0) {
            PRINTF_DBG("ENDing answer_messages_edges. Caught (R_tot)=%d requests, Answered (R_send)=%d.\n",
                       R_tot.size(),R_send.size());

            int status_of_getstatus = 1;
            for (int i = 0; i < R_send.size(); ++i) {
                status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
                while (status_of_getstatus != 0) {
                    PRINTF_DBG("Failing to test a request, its returning %d\n", status_of_getstatus);
                    status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
                }

                if (status_localreq == 0) {
                    PRINTF_DBG("[AM] We are currently waiting for a send to be completed\n");
                    MPI_Wait(&R_send[i], MPI_STATUS_IGNORE);
                    PRINTF_DBG("[AM] Successfully waited the request to be completed\n");
                }
            }
        }

        std::vector<MPI_Request> R_tot;
        std::vector<MPI_Request> R_send;
        //mssleep(DT);
        ++ticks;
    }
    //printf("---answer_request heartbeat---\n");std::cout<<std::flush;
};



template <int DT, int TIMETOL, int BATCH>
void perform_requests(int NNodes,
                      ReferenceContainer REF,
                      unsigned long N,
                      OpenMPHelper &O) {

    std::uniform_int_distribution<int> gen(92412,732532);
    unsigned int SEED = std::stoi(std::getenv("SEED"));
    std::mt19937 rng(SEED);
    int randi = gen(rng);

    long ix;
    int total_processed=0;
    int atomic_helper;

    bool globalstatus = true;
    bool ix_update = false;
    long current_ix;

    bool all_sent_locally = false;
    int tot_locals = 0;
    int sent_locally = 0;

    // Insert here required initialization
    std::list<long> all_indexes_seen;
    if (BATCH<=0)error_report(min_batch_msg);
    bool waiting = false;
    bool bypass = false;
    bool was_last_appearance = false;
    int special_index = -1;
    int rstatus = 0, sstatus = 0;
    MPI_Request requests_send[BATCH] = {MPI_REQUEST_NULL};
    int fsend[BATCH] = {0};
    int frecv[BATCH] = {0};
    int areRecv[BATCH] = {0};
    int myProbes[BATCH] = {0};
    MPI_Request requests_recv[BATCH] = {MPI_REQUEST_NULL};
    InfoVecElem results[BATCH];
    std::queue<int> QAvailable;
    std::list<int> QPend;
    for (int i = 0; i < BATCH; ++i){
        QAvailable.push(i);
    }
    double vval[BATCH], their_vix[BATCH]; // node receptors
    double vval2[BATCH][2], their_vix2[BATCH][2]; // node receptors
    int ixlist[BATCH]; // this is our new helper in order to process various indexes asynchronously (=^-|)-|-/
    int owner[BATCH];
    int retStatus=1;
    int status_rstatus=1;
    int counter=0,MAX_TRIES=3; // #HYPERPARAMS #HYPERPARAMETERS
    bool NONBLOCKING = true;
    bool HANDSHAKE = false;

#pragma omp atomic read
    atomic_helper = *(REF.p_TOT);
    globalstatus = (atomic_helper < NNodes);

    while (globalstatus) {
        // Perform the main task: ask for the required information!
        if (ix_update) {
            auto itBeg = REF.p_ParHelper->data[ix].MissingA.begin();
            auto itEnd = REF.p_ParHelper->data[ix].MissingA.end();
            tot_locals  = 0;
            sent_locally = 0;
            total_processed = 0;
            int localcounter = 0;

            for (auto _it= itBeg; _it != itEnd; ++_it) {
                auto &thread = *_it;
                for (auto it = thread.begin(); it != thread.end(); it++) {
                    ++tot_locals;
                }
            }

            // race-condition unfriendly? Only extensive testing will help us decide :^|
            for (auto _it= itBeg; _it != itEnd; ++_it) {
                auto &thread = *_it;
                for (auto it =  thread.begin();it != thread.end(); it++){
                    ++localcounter;

                    // Retrieve data
                    owner[QAvailable.front()] = std::get<1>(*it);
                    their_vix[QAvailable.front()] = std::get<2>(*it);
                    // Build the result
                    results[QAvailable.front()] = std::make_tuple(0, // placeholder until we get the correct val
                                                                  std::get<0>(*it), // we are inaugurating this indexing model [:<)
                                                                  owner[QAvailable.front()] * N + their_vix[QAvailable.front()]);
                    printf("[PR] About to ask for one node!\n");std::cout<<std::flush;
                    MPI_Ssend(&their_vix[QAvailable.front()],
                              1,
                              MPI_DOUBLE,
                              owner[QAvailable.front()],
                              VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD);
                    printf("[PR] Asked!\n");std::cout<<std::flush;

                    printf("[PR] About to recv!\n");std::cout<<std::flush;
                    MPI_Recv(&vval[QAvailable.front()],
                             1,
                             MPI_DOUBLE,
                             owner[QAvailable.front()],
                             (int) their_vix[QAvailable.front()],
                             MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    printf("[PR] Correctly recieved!\n");std::cout<<std::flush;


                    // Prepare stuff  for the next iteration
                    QPend.push_back(QAvailable.front());
                    QAvailable.pop();

                    if ((QAvailable.empty()) || (localcounter == tot_locals)) {
                        int target_size = 0;
                        if (localcounter == tot_locals){
                            target_size = 0;
                        } else {
                            target_size = BATCH-1;
                        }
                        while (QPend.size() != target_size) {
                            auto i = QPend.begin();
                            while (i != QPend.end()) {
                                std::get<0>(results[*i]) = vval[*i];
                                special_index = ix;
                                (*REF.p_IntHelper)[special_index].ResultsPendProcess.push_back(results[*i]);
                                QAvailable.push(*i);
                                QPend.erase(i++);
                            }
                        }
                    }
                }
            }

//          *******************THE SAME SHOULD BE DONE FOR EDGES PLEASE********************
//          *******************THE SAME SHOULD BE DONE FOR EDGES PLEASE********************
//          *******************THE SAME SHOULD BE DONE FOR EDGES PLEASE********************






//          *******************NOW FOR EDGES************************************************

//            auto itBegEdges = REF.p_ParHelper->data[ix].MissingB.begin();
//            auto itEndEdges = REF.p_ParHelper->data[ix].MissingB.end();
//            tot_locals  = 0;
//            sent_locally = 0;
//            total_processed = 0;
//            localcounter = 0;
//
//            for (auto _it= itBegEdges; _it != itEndEdges; ++_it) {
//                auto &thread = *_it;
//                for (auto it =  thread.begin();it != thread.end(); it++) {
//                    ++tot_locals;
//                }
//            }
//            PRINTF_DBG("(EDGE CASE) Tot locals for index  %d is %d\n", ix, tot_locals);
//            // race-condition unfriendly? Only extensive testing will help us decide [:^|-/<>
//            for (auto _it= itBegEdges; _it != itEndEdges; ++_it) {
//                PRINTF_DBG("(EDGE CASE) About to ask for the process of local elem ix: %d\n",ix);
//                auto &thread = *_it;
//                PRINTF_DBG("(EDGE CASE) SURVIVED 1 \n");
//                for (auto it =  thread.begin();it != thread.end(); it++){
//                    ++localcounter;
//                    PRINTF_DBG("(EDGE CASE) SURVIVED 2 \n");
//
//                    // Retrieve data
//                    owner[QAvailable.front()] = std::get<1>(*it);
//                    their_vix2[QAvailable.front()][0] = std::get<0>(*it);
//                    their_vix2[QAvailable.front()][1] = std::get<2>(*it);
//                    // Build the result
//                    results[QAvailable.front()] = std::make_tuple(0, // placeholder until we get the correct node val
//                                                                  0, // placeholder until we get the correct edge val
//                                                                  // we are inaugurating this indexing model [:<)
//                                                                  owner[QAvailable.front()] * N + their_vix[QAvailable.front()]);
//                    PRINTF_DBG("(EDGE CASE) SURVIVED 3 \n");
//
//                    // Asynchronously start a synchronized sending operation
//                    if (NONBLOCKING) {
//                        retStatus = 1;
//                        requests_send[QAvailable.front()] = MPI_Request();
//                        while (retStatus != 0) {
//                            retStatus = MPI_Isend(&their_vix2[QAvailable.front()][0],//done  NEW STRUCT [0]
//                                                  2, // NEW done 2
//                                                  MPI_DOUBLE,
//                                                  owner[QAvailable.front()],
//                                                  VERTEXVAL_REQUEST_FLAG,
//                                                  MPI_COMM_WORLD,
//                                                  &requests_send[QAvailable.front()]);
//                        }
//                        fsend[QAvailable.front()] = 0;
//                        MPI_Test(&requests_send[QAvailable.front()],
//                                 &fsend[QAvailable.front()],
//                                 MPI_STATUS_IGNORE);
//                    } else {
//                        PRINTF_DBG("Starting blocking ssend");
//                        std::cout << std::flush;
//                        MPI_Ssend(&their_vix2[QAvailable.front()][0], // NEW STRUCT [0]
//                                  2, // NEW 2
//                                  MPI_DOUBLE,
//                                  owner[QAvailable.front()],
//                                  VERTEXVAL_REQUEST_FLAG,
//                                  MPI_COMM_WORLD);
//                        PRINTF_DBG("Ending blocking ssend");
//                        std::cout << std::flush;
//                        fsend[QAvailable.front()] = 1;
//                    }
//
//                    // fsend will monitor if it has arrived
//                    frecv[QAvailable.front()] = 0;
//                    areRecv[QAvailable.front()] = 0;
//                    myProbes[QAvailable.front()] = 0;
//
//
//                    // Prepare stuff  for the next iteration
//                    QPend.push_back(QAvailable.front());
//                    ixlist[QAvailable.front()] = ix;
//                    QAvailable.pop();
//                    PRINTF_DBG("SURVIVED 9 \n");
//
//                    // If QAvailable is empty, devote ourselves to asynchronously waiting for
//                    // answers to arrive. This indirectly means that our QPend has reached length BATCH.
//                    if ((QAvailable.empty()) || (localcounter == tot_locals)) {
//
//                        // Define a target  size that QPend is required to meet.
//                        // a progressive emptying could be coded here :^)
//                        int target_size = 0;
//                        if (localcounter == tot_locals){
//                            // If it is the last lap for this index,
//                            // force the emptying of the container :-)
//                            target_size = 0;
//                        } else {
//                            // If it is not the last lap, just free one
//                            // element.
//                            target_size = BATCH-1;
//                        }
//
//
//                        while (QPend.size() != target_size) {
//
//                            // Iterate through the indexes FIFO style, i.e. std::list travelling
//                            auto i = QPend.begin();
//                            while (i != QPend.end()) {
//
//                                if (fsend[*i] != 1) {
//                                    PRINTF_DBG("Case 1\n");std::cout<<std::flush;
//                                    // If the message has not arrived
//                                    // Update the value to check if it has not arrived
//                                    MPI_Test(&requests_send[*i],
//                                             &fsend[*i],
//                                             MPI_STATUS_IGNORE);
//                                    PRINTF_DBG("Case 1.1\n");std::cout<<std::flush;
//                                    if (fsend[*i] != 1){
//                                        PRINTF_DBG("Case 1.2.1\n");std::cout<<std::flush;
//                                        // If it still has not arrived, we will cancel it and resend it
//                                        retStatus = MPI_Cancel(&requests_send[*i]);
//                                        if (retStatus!=0) PRINTF_DBG("[warning] MPI_Cancel returned %d\n",retStatus);
//                                        MPI_Request_free(&requests_send[*i]);
//                                        requests_send[*i] = MPI_Request();
//                                        PRINTF_DBG("Case 1.2.2\n");std::cout<<std::flush;
//                                        retStatus = MPI_Isend(&their_vix2[*i][0], // NEW STRUCTURE
//                                                              2, // NEW 2
//                                                              MPI_DOUBLE,
//                                                              owner[*i],
//                                                              VERTEXVAL_REQUEST_FLAG,
//                                                              MPI_COMM_WORLD,
//                                                              &requests_send[*i]);
//                                        if (retStatus!=0) PRINTF_DBG("[warning] MPI_Isend returned %d\n",retStatus);
//                                        PRINTF_DBG("Case 1.2.3\n");std::cout<<std::flush;
//                                    }
//                                    PRINTF_DBG("Case 1.3\n");std::cout<<std::flush;
//                                }
//                                if ((fsend[*i] == 1) && (frecv[*i] != 1)){
//                                    PRINTF_DBG("Case 2\n");std::cout<<std::flush;
//                                    // If it has now arrived, check if there is a started reception
//                                    if (areRecv[*i] != 1) {
//                                        PRINTF_DBG("Case 2.1.1\n");std::cout<<std::flush;
//                                        //----------------------------
//                                        myProbes[*i] = 1; // hardcoded
//                                        //----------------------------
//                                        counter = 0;
//                                        while ((myProbes[*i] != 1) && (counter < MAX_TRIES)){
//                                            if (counter == 0) PRINTF_DBG("Case 2.2.1 (first lap)\n");std::cout<<std::flush;
//                                            retStatus = 1;
//                                            while (retStatus != 0){
//                                                retStatus = MPI_Iprobe(owner[*i],
//                                                                       (int) (OFFSET + their_vix2[*i][0]), // done
//                                                                       //(int) their_vix[*i], // NEW TAG
//                                                                       MPI_COMM_WORLD,
//                                                                       &myProbes[*i],
//                                                                       MPI_STATUS_IGNORE);
//                                            }
//                                            counter++;
//                                            //mssleep(DT);
//                                        }
//                                        PRINTF_DBG("Case 2.3\n");std::cout<<std::flush;
//                                        if (myProbes[*i] == 1){
//                                            retStatus = 1;
//                                            PRINTF_DBG("Case 2.3.1\n");std::cout<<std::flush;
//                                            requests_recv[*i] = MPI_Request();
//                                            PRINTF_DBG("Case 2.3.2\n");std::cout<<std::flush;
//                                            recv_nonblocking2(owner[*i], // done NEW FUNCTION
//                                                             requests_recv[*i],
//                                                             vval2[*i][0], // done NEW CONTAINER X2
//                                                              (int) (OFFSET + their_vix2[*i][0])); //NEW TAG: done
//                                                             //(int) their_vix[*i]);
//                                            areRecv[*i] = 1;
//                                            PRINTF_DBG("Case 2.3.3\n");std::cout<<std::flush;
//                                        }
//                                        PRINTF_DBG("Case 2.1-2-3\n");std::cout<<std::flush;
//                                    } else if (areRecv[*i] == 1){
//                                        PRINTF_DBG("Case 2.5\n");std::cout<<std::flush;
//                                        MPI_Test(&requests_recv[*i],
//                                                 &frecv[*i],
//                                                 MPI_STATUS_IGNORE);
//                                        PRINTF_DBG("Case 2.4.1\n");std::cout<<std::flush;
//                                        counter = 0;
//                                        while ((frecv[*i] != 1) && (counter < MAX_TRIES)){
//                                            //PRINTF_DBG("Case 2.4.2\n");std::cout<<std::flush;
//                                            MPI_Test(&requests_recv[*i],
//                                                     &frecv[*i],
//                                                     MPI_STATUS_IGNORE);
//                                            counter++;
//                                            //mssleep(DT);
//                                        }
//                                        PRINTF_DBG("Case 2.5\n");std::cout<<std::flush;
//                                        if (frecv[*i] != 1){
//                                            PRINTF_DBG("Case 2.6.1\n");std::cout<<std::flush;
//                                            // If it still has not arrived, we will cancel it and resend it
//                                            retStatus = MPI_Cancel(&requests_recv[*i]);
//                                            if (retStatus!=0) PRINTF_DBG("[warning] MPI_Cancel returned %d\n",retStatus);
//                                            MPI_Request_free(&requests_recv[*i]);
//                                            frecv[*i] = 0;
//                                            myProbes[*i] = 0;
//                                            PRINTF_DBG("Case 2.6.2\n");std::cout<<std::flush;
//                                        }
//                                    }
//                                }
//
//                                if ((fsend[*i] == 1) && (frecv[*i] == 1)){
//                                    PRINTF_DBG("(EDGES) Case 3\n");std::cout<<std::flush;
//                                    // If it has arrived, store it
//                                    std::get<0>(results[*i]) = vval2[*i][0]; // DONE NEW: NEW X2 STRUCTURE AND RESULTS
//                                    std::get<1>(results[*i]) = vval2[*i][1];
//                                    special_index = ixlist[*i];
//                                    (*REF.p_IntHelper)[special_index].ResultsPendProcess.push_back(results[*i]);
//
//                                    // Debug station // NEW STRUCTURE
//                                    PRINTF_DBG("(EDGES) Recieved a response to our request: val %f %f\n", vval2[*i][0], vval2[*i][1]); std::cout << std::flush;
//
//                                    // Clean
//                                    fsend[*i] = 0;
//                                    frecv[*i] = 0;
//                                    if (NONBLOCKING) MPI_Request_free(&requests_send[*i]);
//                                    MPI_Request_free(&requests_recv[*i]);
//                                    QAvailable.push(*i);
//                                    // Erasing list index
//                                    QPend.erase(i++);
//
//
//
//                                } else { // Index increase for the entire inner loop
//                                    ++i;
//                                }
//
//                                // End of one iteration
//                            }  // End of all iterations
//                            // Now "QPend.size() != target_size" will determine if the while continues.
//                            PRINTF_DBG("About to compute 'QPend.size() != target_size', which yields: %d\n", QPend.size() != target_size); std::cout<< std::flush;
//                            // Optional:
//                            //if (QPend.size() == BATCH) answer_messages_edges<0, 1, BATCH>(REF, O.MY_THREAD_n);
//                        }
//                    }
//                    // Debug station
//                    PRINTF_DBG("DISPATCHED CORRECTLY\n");
//                    std::cout << std::flush;
//                }
//            }
//          *******************END OF FOR EDGES*********************************************


#pragma critical
            {
                REF.p_READY_FOR_INTEGRATION->push(ix);
            }
            ++total_processed;

        } else {
            // *************THIS IS WHAT HAPPENS IF THERE WAS NO INDEX UPDATE****************
            tot_locals  = 0;
            sent_locally = 0;
            total_processed = 0;
            //answer_messages<DT, TIMETOL, BATCH>(REF, O.MY_THREAD_n); // #HYPERPARAMS #HYPERPARAMETERS
        }
        ix_update = false;
#pragma omp critical
        {
            if (!REF.p_CHECKED->empty()) {
                ix = REF.p_CHECKED->front();
                REF.p_CHECKED->pop();
                ix_update = true;
            }
        }
        // Add to the number of totals processed the ones from previous lap :-)
        if (total_processed != 0) {
#pragma omp atomic update
            *(REF.p_TOT) += total_processed;
        }

        // if there was no available index, check if it is still worth looping.
        if (!ix_update) {
#pragma omp atomic read
            atomic_helper = *(REF.p_TOT);
            globalstatus = (atomic_helper < NNodes); // this should be a strict equality but its buggy :^(
        }
        PRINTF_DBG("TOT=%d, NNodes=%d, ix_update=%d, total_processed=%d, ix=%d\n",
                   atomic_helper, NNodes, ix_update, total_processed, ix);std::cout<<std::flush;
    }
    PRINTF_DBG("Final termination of perform_requests :-)\n");std::cout<<std::flush;

} // end of function




void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H