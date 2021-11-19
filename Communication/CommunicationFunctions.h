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
void build_answer_edges(double * answer, ReferenceContainer &REF, double * ix, int owner, int MyNProc);

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
    MPI_Message M;
    MPI_Status  S;
    MPI_Request R;
    int TRIES=0;
    int t=0;
    bool firstlap=true;
    while (TRIES < 3){
        MPI_Status status;
        if ((flag == 1) || firstlap) {
            flag = 0;
            buffer = -9995.0;
            R = MPI_Request();
            M = MPI_Message();
            S = MPI_Status();

            MPI_Improbe(MPI_ANY_SOURCE, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD,
                    &flag, &M, &status);
        }
        if (firstlap) firstlap = false;

        while ((flag != 1) && (t<TIMETOL)){

            MPI_Improbe(MPI_ANY_SOURCE, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD,
                        &flag, &M, &status);
            ++t;
            mssleep(DT);
        }
        if (t>=TIMETOL) ++TRIES;
        t=0;
        if (flag == 1){
            MPI_Mrecv(&buffer, 1, MPI_DOUBLE, &M, &S);
            build_answer(answer, REF, buffer, S.MPI_SOURCE, MYPROC);
            PRINTF_DBG("About to send one!\n");std::cout<<std::flush;
            MPI_Ssend(&answer, 1, MPI_DOUBLE, S.MPI_SOURCE, (int) buffer, MPI_COMM_WORLD);
            PRINTF_DBG("Correctly sent one!\n");std::cout<<std::flush;
        }
    }
};






template<int DT, int TIMETOL, int BATCH>
void answer_messages_edges(ReferenceContainer &REF,int MYTHR) {

    int MYPROC = REF.p_ComHelper->WORLD_RANK[MYTHR];
    int NPROCS = REF.p_ComHelper->WORLD_SIZE[MYTHR];
    int NDISPATCHED = 0;
    int flag;
    int probe_status = 1;
    double ix;
    int statusreq_status = 1;
    double buffer[2];
    double answer[2];
    int ticks = 0;
    int statusFree = 0;
    int NTOT = 0;
    int NERR = 0;
    MPI_Message M;
    MPI_Status  S;
    MPI_Request R;
    int TRIES=0;
    int t=0;
    bool firstlap=true;
    while (TRIES < 3){
        MPI_Status status;
        if ((flag == 1) || firstlap) {
            flag = 0;
            R = MPI_Request();
            M = MPI_Message();
            S = MPI_Status();

            MPI_Improbe(MPI_ANY_SOURCE, EDGEVAL_REQUEST_FLAG, MPI_COMM_WORLD,
                        &flag, &M, &status);
        }
        if (firstlap) firstlap = false;

        while ((flag != 1) && (t<TIMETOL)){

            MPI_Improbe(MPI_ANY_SOURCE, EDGEVAL_REQUEST_FLAG, MPI_COMM_WORLD,
                        &flag, &M, &status);
            ++t;
            mssleep(DT);
        }
        if (t>=TIMETOL) ++TRIES;
        t=0;
        if (flag == 1){
            MPI_Mrecv(&buffer[0], 2, MPI_DOUBLE, &M, &S);
            build_answer_edges(&answer[0], REF, &buffer[0], S.MPI_SOURCE, MYPROC);
            PRINTF_DBG("About to send one!\n");std::cout<<std::flush;
            MPI_Ssend(&answer[0],
                      2,
                      MPI_DOUBLE,
                      S.MPI_SOURCE,
                      (int) ((int) buffer[1] + (int) OFFSET),
                      MPI_COMM_WORLD);
            PRINTF_DBG("Correctly sent one!\n");std::cout<<std::flush;
        }
    }
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
    unsigned long uix;
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
                    their_vix[QAvailable.front()] = (double) std::get<2>(*it);
                    // Build the result
                    results[QAvailable.front()] = std::make_tuple((double) 0, // placeholder until we get the correct val
                                                                  (double) std::get<0>(*it), // we are inaugurating this indexing model [:<)
                                                      owner[QAvailable.front()] * N + (int) their_vix[QAvailable.front()]);
                    PRINTF_DBG("[PR] About to ask for one node!\n");std::cout<<std::flush;
                    MPI_Ssend(&their_vix[QAvailable.front()],
                              1,
                              MPI_DOUBLE,
                              owner[QAvailable.front()],
                              VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD);
                    PRINTF_DBG("[PR] Asked!\n");std::cout<<std::flush;

                    PRINTF_DBG("[PR] About to recv!\n");std::cout<<std::flush;
                    MPI_Recv(&vval[QAvailable.front()],
                             1,
                             MPI_DOUBLE,
                             owner[QAvailable.front()],
                             (int) their_vix[QAvailable.front()],
                             MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    PRINTF_DBG("[PR] Correctly recieved!\n");std::cout<<std::flush;


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

            auto itBegE = REF.p_ParHelper->data[ix].MissingB.begin();
            auto itEndE = REF.p_ParHelper->data[ix].MissingB.end();
            tot_locals  = 0;
            sent_locally = 0;
            total_processed = 0;
            localcounter = 0;

            for (auto _it= itBegE; _it != itEndE; ++_it) {
                auto &thread = *_it;
                for (auto it = thread.begin(); it != thread.end(); it++) {
                    ++tot_locals;
                }
            }



            // race-condition unfriendly? Only extensive testing will help us decide :^|
            for (auto _it= itBegE; _it != itEndE; ++_it) {
                auto &thread = *_it;
                for (auto it =  thread.begin();it != thread.end(); it++){
                    ++localcounter;


                    // Retrieve data
                    owner[QAvailable.front()] = (int) std::get<1>(*it);
                    their_vix2[QAvailable.front()][0] = (double) std::get<0>(*it);
                    their_vix2[QAvailable.front()][1] = (double) std::get<2>(*it);
                    // Build the result
                    results[QAvailable.front()] = std::make_tuple(0.0, // placeholder until we get the correct node val
                                                                  0.0, // placeholder until we get the correct edge val
                            // we are inaugurating this indexing model [:<)
                                                                  owner[QAvailable.front()] * N + (int) their_vix2[QAvailable.front()][1]);


                    PRINTF_DBG("[PR] About to ask for one node and edge!\n");std::cout<<std::flush;
                    PRINTF_DBG("About to send ixs %f and %f\n",their_vix2[QAvailable.front()][0], their_vix2[QAvailable.front()][1]);std::cout<<std::flush;
                    MPI_Ssend(&their_vix2[QAvailable.front()][0],
                              2,
                              MPI_DOUBLE,
                              owner[QAvailable.front()],
                              EDGEVAL_REQUEST_FLAG, MPI_COMM_WORLD);
                    PRINTF_DBG("[PR] Asked!\n");std::cout<<std::flush;

                    PRINTF_DBG("[PR] About to recv!\n");std::cout<<std::flush;
                    MPI_Recv(&vval2[QAvailable.front()][0],
                             2,
                             MPI_DOUBLE,
                             owner[QAvailable.front()],
                             (int) (int) (OFFSET + (int) their_vix2[QAvailable.front()][1]),
                             MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    PRINTF_DBG("[PR] Correctly recieved!\n");std::cout<<std::flush;

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
                                std::get<0>(results[*i]) = vval2[*i][0];
                                std::get<1>(results[*i]) = vval2[*i][1];
                                special_index = ix;
                                (*REF.p_IntHelper)[special_index].ResultsPendProcess.push_back(results[*i]);
                                QAvailable.push(*i);
                                QPend.erase(i++);
                            }

                        }
                    }
                }
            }


#pragma critical
            {
                (*REF.p_READY_FOR_INTEGRATION).first.push(ix);
                (*REF.p_READY_FOR_INTEGRATION).second.push(uix);
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
            if (!(*REF.p_CHECKED).first.empty()) { // Assuming they have both the same length
                ix = (*REF.p_CHECKED).first.front();
                uix = (*REF.p_CHECKED).second.front();
                (*REF.p_CHECKED).first.pop();
                (*REF.p_CHECKED).second.pop();
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
        if (globalstatus && (!ix_update)) mssleep(5);
    }
    PRINTF_DBG("Final termination of perform_requests :-)\n");std::cout<<std::flush;

} // end of function




void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H