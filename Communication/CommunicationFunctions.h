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

#ifndef VERTEXVAL_REQUEST_FLAG
#define VERTEXVAL_REQUEST_FLAG -1
#endif

int destroyRequestReturnInteger(MPI_Request &R);
void destroyRequest(MPI_Request &R, int &NERR);
void destroyRequestWithoutCounter(MPI_Request &R);
void send_nonblocking(int owner, MPI_Request &r, double &ix, int TAG);
void recv_nonblocking(int owner, MPI_Request &r, double &result, int TAG);

void irespond_value(ReferenceContainer &REF, double ix, int owner, MPI_Request & R, int MyNProc);




template<int DT, int TIMETOL, int BATCH>
void answer_messages(ReferenceContainer &REF,int MYTHR){

    int MYPROC = REF.p_ComHelper->WORLD_RANK[MYTHR];
    std::vector<MPI_Request> R_tot;
    std::vector<MPI_Request> R_send;
    R_tot.push_back(MPI_Request());
    R_send.push_back(MPI_Request());
    int NDISPATCHED = 0;

    // lay the probes for all the p rocs
    int flagprobes[MYTHR];
    for (int i=0; i<MYTHR; i++){
        flagprobes[i] = 0;
        if (i != MYPROC){
            MPI_Iprobe(i,VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
        }
    }

    // initialize auxiliary variables
    int ticks = 0;
    std::set<int> answered;
    answered.insert(MYPROC);
    double ix;
    int i = 0;
    int NTOT=0;
    int status = 0;
    int status_localreq = 0;
    int NERR = 0;

    for (int tk=0; tk<TIMETOL; ++tk){
        //printf("running outer loop\n");
        for (int j=0; j < BATCH; ++j) {
            //printf("running inner loop\n");
            if (flagprobes[i] == 1) { // a message is available! recieve it!
                PRINTF_DBG("Entered flagprobes\n");
                MPI_Irecv(&ix,
                          1,//count
                          MPI_DOUBLE, // type
                          i, // destination
                          VERTEXVAL_REQUEST_FLAG, // this is the flag, meaning a request for a vertex value :-)
                          MPI_COMM_WORLD, &R_tot[R_tot.size()-1]);
                MPI_Request_get_status(R_tot[R_tot.size()-1], &status_localreq, MPI_STATUS_IGNORE);
                if (status_localreq==1){ // If we were first to capture the message, proceed.
                    PRINTF_DBG("We effectively captured a vertex info request :-)\n");
                    irespond_value(REF, ix, i, R_send[R_send.size()-1], MYPROC);
                    ++NTOT;
                    PRINTF_DBG("We effectively answered asynchronously a vertex info request :-)\n");
                    R_send.push_back(MPI_Request());
                    NDISPATCHED++;
                } else {
                    printf("We were faced with a probe that indicated an incoming message but we couldnt capture it :o\n");
                }
                status_localreq = 0;
                flagprobes[i] = 0;
                answered.insert(i); // it was answered, just not by us :-)
                //destroyRequest(R_tot[R_tot.size()-1], NERR);
                R_tot.push_back(MPI_Request());
            }
            if (answered.size() == MYTHR){
                // we should just refresh answered :-)
                answered = std::set<int>();
                answered.insert(MYPROC);
            }
            PRINTF_DBG("Arrived to C\n");
            ++i;
            if (i >= MYTHR) {
                PRINTF_DBG("Entered D\n");
                i = 0;
            }
            while (answered.count(i) == 1){
                //printf("Entered E\n");
                ++i;
                if (i >= MYTHR) {
                    i = 0;
                }
            }

        }
        mssleep(DT);
        ++ticks;
    }
    //printf("exited :-)");
    for (int i=0; i<R_tot.size()-1;++i){
        PRINTF_DBG("Currently at GA\n");
        //destroyRequest(R_tot[i], NERR);
    }
    for (int i=0; i<R_send.size()-1;++i){ // we wait for our answers to arrive before descoping
        MPI_Request_get_status(R_send[i], &status_localreq, MPI_STATUS_IGNORE);
        if (status_localreq == 0){
            printf("We are currently waiting for a send to be completed :-(\n");
            std::cout << std::flush;
            MPI_Wait(&R_send[i], MPI_STATUS_IGNORE);
        }
        //destroyRequest(R_send[i], NERR);
        PRINTF_DBG("Successfully answered one request! ^.^\n");
    }
    PRINTF_DBG("\nExiting AnswerMsg: there were %d failures to erase requests.\n", NERR);
    //printf("Exiting answer_messages, there were %d responded messages\n",NTOT);
};




template <int BATCH>
void GetOneMsg(int ix,
               ReferenceContainer &REF,
               unsigned long N){

    if (BATCH<=0)error_report(min_batch_msg); // guard: we need a batch size of at least 1 :-)
    printf("entered getonemsg...\n");

    bool waiting = false;
    int rstatus = 0, sstatus = 0;
    MPI_Request requests_send[BATCH], requests_recv[BATCH];
    InfoVecElem results[BATCH];
    std::queue<int> QAvailable;
    std::list<int> QPend;
    for (int i = 0; i < BATCH; ++i){
        QAvailable.push(i);
    }

    double vval[BATCH], their_vix[BATCH];
    int owner[BATCH];

    // Collect missing nodes:
    // Access to "MissingA" can be direct as there are not
    // race conditions :-)
    auto itBeg = REF.p_ParHelper->data[ix].MissingA.begin();
    auto itEnd = REF.p_ParHelper->data[ix].MissingA.end();
    for (auto _it= itBeg; _it != itEnd; ++_it) { // race-condition unfriendly?
        auto &thread = *_it;
        printf("There were indeed  elements to iterate :-)\n");
        std::cout << "Length of this section is: " << itEnd - itBeg << std::endl;
        std::cout << std::flush; // DEBUGGING
        for (auto it =  thread.begin();it != thread.end(); it++){
            // retrieve data
            owner[QAvailable.front()] = std::get<1>(*it);
            their_vix[QAvailable.front()] = std::get<2>(*it);
            // send request for missing info
            printf("Sending and recieving (nonblocking): asking %d for vertex w index: %d \n",
                   owner[QAvailable.front()], their_vix[QAvailable.front()]);
            std::cout  << std::flush; // DEBUGGING
            send_nonblocking(owner[QAvailable.front()], // send a request (tag=-1) for some ix's vertex value
                                     requests_send[QAvailable.front()],
                                     their_vix[QAvailable.front()], VERTEXVAL_REQUEST_FLAG);
            recv_nonblocking(owner[QAvailable.front()], // recieve the respone (tag=ix) for that  ix
                                     requests_recv[QAvailable.front()], // flag = (int) index :-)
                                     vval[QAvailable.front()], (int) their_vix[QAvailable.front()]);
            results[QAvailable.front()] = std::make_tuple(0, // placeholder until we get the correct val
                                                          std::get<0>(*it),
                                                          owner[QAvailable.front()] * N + their_vix[QAvailable.front()]);
            QPend.push_back(QAvailable.front());
            QAvailable.pop();
            if (QAvailable.empty()) {
                waiting = true;
                while (waiting) {
                    printf("We are waiting :-(\n");
                    std::cout << std::flush;
                    mssleep(100);
                    auto i = QPend.begin();
                    while (i != QPend.end()) {
                        MPI_Request_get_status(requests_send[*i], &sstatus, MPI_STATUS_IGNORE);
                        MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                        if ((sstatus==1) && (rstatus==1)) {
                            printf("Recieved a response to our request: val %f\n", vval[*i]);
                            std::cout << std::flush;
                            std::get<0>(results[*i]) = vval[*i];
                            // No lock is needed: this index is entirely ours :-)
                            (*REF.p_IntHelper)[ix].ResultsPendProcess.push_back(results[*i]);
                            QAvailable.push(*i);
                            QPend.erase(i++);
                            waiting = false;
                            requests_send[*i] = MPI_Request();
                            requests_recv[*i] = MPI_Request();
                        } else {
                            ++i;
                        }
                    }
                }
            }
        }
    }
    // TIME FOR EDGES NOW...
    // Collect missing nodes+edges
    //****************************
    //The following section is commented
    // as to debug only the previous one :-)
    //****************************
//    for (auto & thread : P.data[ix].MissingB) {
//        for (auto it =  thread.begin();
//             it != thread.end(); it++){
//            // retrieve data
//            double vval, eval;
//            int owner = std::get<1>(*it);
//            int their_vix = std::get<2>(*it);
//            // ask for it, recieve it, and store it.
//            //ask_for_node_and_vertex(owner, vval, eval, H, their_vix, g);
//#pragma omp atomic
//            I[ix].ResultsPendProcess.emplace_back(vval,eval,owner * N + their_vix);
//        }
//    }

// BEFORE EXITING PLEASE process the remaining cases :-) N :-)
    waiting = (QPend.size()>0);
    while (waiting) {
        printf("We are waiting before exiting, i.e. CLEANING! :-(\n");
        std::cout << std::flush;
        mssleep(100);
        auto i = QPend.begin();
        while (i != QPend.end()) {
            MPI_Request_get_status(requests_send[*i], &sstatus, MPI_STATUS_IGNORE);
            MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
            if ((sstatus==1) && (rstatus==1)) {
                printf("Recieved a response to our request: val %f\n", vval[*i]);
                std::cout << std::flush;
                std::get<0>(results[*i]) = vval[*i];
                // No lock is needed: this index is entirely ours :-)
                (*REF.p_IntHelper)[ix].ResultsPendProcess.push_back(results[*i]);
                QAvailable.push(*i);
                QPend.erase(i++);
                waiting = false;
                requests_send[*i] = MPI_Request();
                requests_recv[*i] = MPI_Request();
            } else {
                ++i;
            }
        }
    }
printf("Exiting get one msg!\n");
}






template <int DT, int TIMETOL, int BATCH>
void perform_requests(int NNodes,
                       ReferenceContainer REF,
                       unsigned long N,
                       OpenMPHelper &O) {

    long ix;
    int total_processed=0;
    int atomic_helper;

    bool globalstatus = true;
    bool ix_update = false;

#pragma omp atomic read
    atomic_helper = *(REF.p_TOT);
    globalstatus = (atomic_helper < NNodes);

    while (globalstatus) {
        // Perform the main task: ask for the required information!
            if (ix_update) {
                //mssleep(1000);
                //std::cout << "[DEBUG] about to call GetOneMsg" << std::endl;
                printf("About to ask for the process of local elem ix: %d\n",ix);
                GetOneMsg<BATCH>(ix, REF, N);
                PRINTF_DBG("Skipped sending messages :-)\n");
#pragma omp atomic update
                ++ *(REF.p_TOT);
                ++total_processed;
                PRINTF_DBG("increasing TOT, proc \n");
                // Also communicate to other threads that there is one new item pending integration :o
#pragma critical
                {
                    REF.p_READY_FOR_INTEGRATION->push(ix);
                }
                printf("Reporting so far having requested %d nodes... there are in total %d nodes\n",
                       total_processed, atomic_helper);
                ix_update = false;
            } else {
                answer_messages<DT, TIMETOL, BATCH>(REF, O.MY_THREAD_n);
            }

#pragma omp critical
{
            if (!REF.p_CHECKED->empty()) {
                ix = REF.p_CHECKED->front();
                REF.p_CHECKED->pop();
                ix_update = true;
            }
}

#pragma omp atomic read // check if this loop still makes sense
            atomic_helper = *(REF.p_TOT);
            globalstatus = (atomic_helper < NNodes);
        }
} // end of function















void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
