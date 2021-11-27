//
// Created by m4zz31 on 5/11/21.
//

#ifndef CPPPROJCT_COMMUNICATIONFUNCTIONS_H
#define CPPPROJCT_COMMUNICATIONFUNCTIONS_H

#include <set>
#include <iterator>
#include <random>

#include "../macros/macros.h"

#include "../GraphClasses/GeneralGraph.h"
#include "../GraphClasses/GraphFunctions.h"

#include "../Utils/HelperClasses.h"
#include "../Utils/error.h"
#include "../Utils/global_standard_messages.h"
#include "../Utils/msleep.h"



void build_answer_nodes(double &answer, ReferenceContainer &REF, double ix, int owner, int &MyNProc);
void build_answer_edges(double * answer, ReferenceContainer &REF, double * ix, int owner, int &MyNProc);

template <int T>
void build_answer_generic(double * answer, ReferenceContainer &REF, double * ix, int owner, int &MyNProc){
    if (T == 0){
        build_answer_nodes(*answer, REF, *ix, owner, MyNProc);
    } else if (T == 1) {
        build_answer_edges(answer, REF, ix, owner, MyNProc);
    }
};



// Trying to dispatch requests for Runge Kutta terms :-)
template<int DT, int TIMETOL, int BATCH>
void answer_field_requests(ReferenceContainer &REF,int MYTHR, int fieldOrder){
    // Gonna receive requests for runge-kutta terms :-)
    // REF can access the RK4 terms via (*REF.p_RK4)[ix]
    // We should get requests for specific terms like
    // TAG for RK1 request, TAG for RK1 answer, etc...
    // the mssg says: [our ix] (int) and we respond with [our ix, rk value :-)]
    int flag[4] = {0};
    int TRIES=0;
    MPI_Request R;
    MPI_Message M;
    MPI_Status S;
    int buffer[2];
    double answer;
    int t=0;
    bool firstlap[4]= {true};
    int ASKING_TAGS[4] = {K1_REQUEST, K2_REQUEST, K3_REQUEST, K4_REQUEST};
    int ANSWERING_TAGS[4] = {K1_ANSWER, K2_ANSWER, K3_ANSWER, K4_ANSWER};
    while (TRIES < 4){
        for (int i = fieldOrder-1-1; i < fieldOrder-1; ++i) {
            flag[i] = 0;
            MPI_Status status;
            if ((flag[i] == 1) || firstlap[i]) {
                flag[i] = 0;
                buffer[0] = -9995;
                buffer[1] = -9994;
                R = MPI_Request();
                M = MPI_Message();
                S = MPI_Status();
                MPI_Improbe(MPI_ANY_SOURCE, ASKING_TAGS[i], MPI_COMM_WORLD,
                            &flag[i], &M, &status);
            }
            if (firstlap[i]) firstlap[i] = false;
            while ((flag[i] != 1) && (t < TIMETOL)) {

                MPI_Improbe(MPI_ANY_SOURCE, ASKING_TAGS[i], MPI_COMM_WORLD,
                            &flag[i], &M, &status);
                ++t;
                mssleep(DT);
            }
            if (t >= TIMETOL) ++TRIES;
            t = 0;
            if (flag[i] == 1) {
                MPI_Mrecv(&buffer, 2, MPI_INT, &M, &S);
                PRINTF_DBG("[afr] WE JUST RECEIVED :-) buffer is %d and %d\n", buffer[0], buffer[1]);
                bool isReady = false;
                if (i==0){
#pragma omp atomic read
                    isReady = REF.p_LayHelper->data[buffer[1]].RK1_status;
                    while (!isReady){
                        PRINTF_DBG("RK1 status not found, keeping going...");std::cout<<std::flush;
                        mssleep(DT);
#pragma omp atomic read
                        isReady = REF.p_LayHelper->data[buffer[1]].RK1_status;
                    }
#pragma omp atomic read
                    answer = REF.p_LayHelper->data[buffer[1]].RK1[0];
                } else if (i==1) {
#pragma omp atomic read
                    isReady = REF.p_LayHelper->data[buffer[1]].RK2_status;
                    while (!isReady){
                        PRINTF_DBG("RK2 status not found, keeping going...");std::cout<<std::flush;
                        mssleep(DT);
#pragma omp atomic read
                        isReady = REF.p_LayHelper->data[buffer[1]].RK2_status;
                    }
#pragma omp atomic read
                    answer = REF.p_LayHelper->data[buffer[1]].RK2[0];
                } else if (i==2) {
#pragma omp atomic read
                    isReady = REF.p_LayHelper->data[buffer[1]].RK3_status;
                    while (!isReady){
                        PRINTF_DBG("RK3 status not found, keeping going...");std::cout<<std::flush;
                        mssleep(DT);
#pragma omp atomic read
                        isReady = REF.p_LayHelper->data[buffer[1]].RK3_status;
                    }

#pragma omp atomic read
                    answer = REF.p_LayHelper->data[buffer[1]].RK3[0];
                } else if (i==3) {
#pragma omp atomic read
                    isReady = REF.p_LayHelper->data[buffer[1]].RK4_status;
                    while (!isReady){
                        PRINTF_DBG("RK4 status not found, keeping going...");std::cout<<std::flush;
                        mssleep(DT);
#pragma omp atomic read
                        isReady = REF.p_LayHelper->data[buffer[1]].RK4_status;
                    }

#pragma omp atomic read
                    answer = REF.p_LayHelper->data[buffer[1]].RK4[0];
                }
                PRINTF_DBG("[afr] About to send :-)! tag %d\n", buffer[0]);std::cout<<std::flush;
                MPI_Ssend(&answer, 1, MPI_DOUBLE, S.MPI_SOURCE, buffer[0], MPI_COMM_WORLD);
                PRINTF_DBG("[afr] Finished sending :-)! tag %d\n", buffer[0]);
            }
        }
    }
}





template <typename SpecificRequestObject>
class RequestObject : public SpecificRequestObject {
    public:
        RequestObject(int type): SpecificRequestObject(type){};
        RequestObject(){};
};

class BaseSpecificRequestObject{
public:
    int recvLength;
    int sendLength;
    int recvTag;
    int sendTag;
};

class FieldRequestObject : public BaseSpecificRequestObject{
public:
    int field=-1;
    int recvLength = 2;
    int  sendLength = 1;
    FieldRequestObject(int fieldOrder): field(fieldOrder){};
    void buildSendTag(int * data);
    void buildRecvTag(int * data);
    void computeAnswer(ReferenceContainer &REF, int * buffer, double * answer);
    void computeReady(ReferenceContainer &REF, int * buffer, bool &isReady);
};



template<int DT, int TIMETOL, int BATCH, typename RequestClass>
void generic_answer_requests(ReferenceContainer &REF,int MYTHR, RequestClass ReqObj){
    int flag = 0;
    int TRIES=0;
    MPI_Request R;
    MPI_Message M;
    MPI_Status S;
    int buffer[ReqObj.recvLength];
    double answer[ReqObj.sendLength];
    int t=0;
    bool firstlap = true;
    ReqObj.buildRecvTag(&buffer[0]);

    while (TRIES < 4){
        flag = 0;
        MPI_Status status;
        if ((flag == 1) || firstlap) {
            flag = 0;
            R = MPI_Request();
            M = MPI_Message();
            S = MPI_Status();
            MPI_Improbe(MPI_ANY_SOURCE,
                        ReqObj.recvTag,
                        MPI_COMM_WORLD,
                        &flag,
                        &M,
                        &status);
        }
        if (firstlap) firstlap = false;
        while ((flag != 1) && (t < TIMETOL)) {
            MPI_Improbe(MPI_ANY_SOURCE,
                        ReqObj.recvTag,
                        MPI_COMM_WORLD,
                        &flag,
                        &M,
                        &status);
            ++t;
            mssleep(DT);
        }
        if (t >= TIMETOL) ++TRIES;
        t = 0;
        if (flag == 1) {
            MPI_Mrecv(&buffer, ReqObj.recvLength, MPI_INT, &M, &S);
            bool isReady = false;

#pragma omp critical
{
            ReqObj.computeReady(REF, buffer, isReady);
};
            while (!isReady) {
#pragma omp critical
{
                ReqObj.computeReady(REF, buffer, isReady);
};
                mssleep(DT);
            }
            ReqObj.computeAnswer(REF, buffer, answer);
            ReqObj.buildSendTag(&buffer[0]);
            MPI_Ssend(&answer[0],
                      ReqObj.sendLength,
                      MPI_DOUBLE,
                      S.MPI_SOURCE,
                      ReqObj.sendTag,
                      MPI_COMM_WORLD);
        }
    }
}















template<int DT, int TIMETOL, int BATCH>
void perform_field_requests(ReferenceContainer &REF,int MYPROC, int fieldOrder,std::queue<long> * queue){

    int ASKING_TAGS[4] = {K1_REQUEST, K2_REQUEST, K3_REQUEST, K4_REQUEST};
    bool keep_working = true;
    long ix;
    int L;
#pragma omp critical
{
    if (!queue->empty()){
        ix = queue->front();
        queue->pop();
    } else {
        keep_working = false;
    }
}
    while (keep_working){
        if (fieldOrder==2){
            L = REF.p_LayHelper->data[ix].RK1.size();
        } else if (fieldOrder==3) {
            L = REF.p_LayHelper->data[ix].RK2.size();
        } else if (fieldOrder==4) {
            L = REF.p_LayHelper->data[ix].RK2.size();
        } else {
            printf("[FATAL] field order requested to perform_field_requests does not exist!\n");std::cout<<std::flush;
            exit(1);
        }

        for (int i=0; i<L-1 ; ++i){
            double recvBuffer;
            if (((int) std::get<2>((*REF.p_IntHelper)[ix].ixMap[i])) == MYPROC) {
                if (fieldOrder == 2) {
                    bool isReadyYet = false; //RECENT: changed from get<2> to get<1>
                    unsigned long owner = std::get<1>((*REF.p_IntHelper)[ix].ixMap[i]);
#pragma omp atomic read
                    isReadyYet = REF.p_LayHelper->data[owner].RK1_status;
                    while (!isReadyYet) {
                        mssleep(DT);
#pragma omp atomic read
                        isReadyYet = REF.p_LayHelper->data[owner].RK1_status;
                    }
#pragma omp atomic read
                    recvBuffer = REF.p_LayHelper->data[owner].RK1[0];
                } else if (fieldOrder == 3) {
                    bool isReadyYet = false;
                    unsigned long owner = std::get<2>((*REF.p_IntHelper)[ix].ixMap[i]);
#pragma omp atomic read
                    isReadyYet = REF.p_LayHelper->data[owner].RK2_status;
                    while (!isReadyYet) {
                        mssleep(DT);
#pragma omp atomic read
                        isReadyYet = REF.p_LayHelper->data[owner].RK2_status;
                    }
#pragma omp atomic read
                    recvBuffer = REF.p_LayHelper->data[owner].RK2[0];
                } else if (fieldOrder == 4) {
                    bool isReadyYet = false;
                    unsigned long owner = std::get<2>((*REF.p_IntHelper)[ix].ixMap[i]);
#pragma omp atomic read
                    isReadyYet = REF.p_LayHelper->data[owner].RK3_status;
                    while (!isReadyYet) {
                        mssleep(DT);
#pragma omp atomic read
                        isReadyYet = REF.p_LayHelper->data[owner].RK3_status;
                    }
#pragma omp atomic read
                    recvBuffer = REF.p_LayHelper->data[owner].RK3[0];
                };
            } else {
                int sendBuffer[2] = {
                        (int) std::get<0>((*REF.p_IntHelper)[ix].ixMap[i]),
                        (int) std::get<1>((*REF.p_IntHelper)[ix].ixMap[i])
                };

                printf("[pfr] About to send! sendbuffer says %d and %d... asking tag is %d\n", sendBuffer[0], sendBuffer[1], ASKING_TAGS[fieldOrder-2]);
                //std::cout << "Originally it is: " << std::get<0>((*REF.p_IntHelper)[ix].ixMap[i]) << std::endl;
                //std::cout << "L is " << L << " and [ix].ixMap's size is: " << (*REF.p_IntHelper)[ix].ixMap.size() << std::endl;
                MPI_Ssend(&sendBuffer,
                          2,
                          MPI_INT,
                          (int) std::get<2>((*REF.p_IntHelper)[ix].ixMap[i]),
                          ASKING_TAGS[fieldOrder-2],
                          MPI_COMM_WORLD);
                PRINTF_DBG("[pfr] About to receive! tag is %d\n",
                       (int) std::get<0>((*REF.p_IntHelper)[ix].ixMap[i])
                );
                MPI_Recv(&recvBuffer,
                         1,
                         MPI_DOUBLE,
                         (int) std::get<2>((*REF.p_IntHelper)[ix].ixMap[i]),
                         (int) std::get<0>((*REF.p_IntHelper)[ix].ixMap[i]),
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                PRINTF_DBG("[pfr] recieved %f!\n",recvBuffer);
            }

            if (fieldOrder == 2) {
                REF.p_LayHelper->data[ix].RK1[i+1] = recvBuffer;
            } else if (fieldOrder == 3) {
                REF.p_LayHelper->data[ix].RK2[i+1] = recvBuffer;
            } else if (fieldOrder == 4) {
                REF.p_LayHelper->data[ix].RK3[i+1] = recvBuffer;
            };
        }
#pragma omp critical
 {
         if (!queue->empty()){
             ix = queue->front();
             queue->pop();
         } else {
             keep_working = false;
         }
 }
     }
};

template<int DT, int TIMETOL, int BATCH>
void answer_messages(ReferenceContainer &REF,int MYTHR) {

    int flag;
    double ix;
    double buffer;
    double answer;
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
            build_answer_generic<0>(&answer, REF, &buffer, S.MPI_SOURCE, REF.p_ComHelper->WORLD_RANK[MYTHR]);
            PRINTF_DBG("About to send one!\n");std::cout<<std::flush;
            MPI_Ssend(&answer, 1, MPI_DOUBLE, S.MPI_SOURCE, (int) buffer, MPI_COMM_WORLD);
            PRINTF_DBG("Correctly sent one!\n");std::cout<<std::flush;
        }
    }
};






template<int DT, int TIMETOL, int BATCH>
void answer_messages_edges(ReferenceContainer &REF,int MYTHR) {

    int flag;
    double ix;
    double buffer[2];
    double answer[2];
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
            build_answer_generic<1>(&answer[0], REF, &buffer[0], S.MPI_SOURCE, REF.p_ComHelper->WORLD_RANK[MYTHR]);
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
    double vval[BATCH], their_vix[BATCH];
    double vval2[BATCH][2], their_vix2[BATCH][2];
    int ixlist[BATCH]; 
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

            for (auto _it= itBeg; _it != itEnd; ++_it) {
                auto &thread = *_it;
                auto it = thread.begin();
                while (it !=  thread.end()) {
                    ++localcounter;

                    owner[QAvailable.front()] = std::get<1>(*it);
                    their_vix[QAvailable.front()] = (double) std::get<2>(*it);
                    results[QAvailable.front()] = std::make_tuple((double) 0, // placeholder until we get the correct val
                                                                  (double) std::get<0>(*it),
                          (unsigned long) (((unsigned long) owner[QAvailable.front()]) *  N + (unsigned long) their_vix[QAvailable.front()] ),
                                                                  (unsigned long) std::get<2>(*it),
                                                                  (unsigned long) std::get<1>(*it));
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
                    thread.erase(it++);
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



            for (auto _it= itBegE; _it != itEndE; ++_it) {
                auto &thread = *_it;
                auto it = thread.begin();
                while (it !=  thread.end()) {
                    ++localcounter;

                    // Retrieve data
                    owner[QAvailable.front()] = (int) std::get<1>(*it);
                    their_vix2[QAvailable.front()][0] = (double) std::get<0>(*it);
                    their_vix2[QAvailable.front()][1] = (double) std::get<2>(*it);
                    // Build the result
                    results[QAvailable.front()] = std::make_tuple(0.0, // placeholder until we get the correct node val
                                                                  0.0, // placeholder until we get the correct edge val
                          (unsigned long) (((unsigned long) owner[QAvailable.front()]) *  N + (unsigned long) their_vix2[QAvailable.front()][1] ),
                           (unsigned long) std::get<2>(*it),
                            (unsigned long) std::get<1>(*it));


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
                    thread.erase(it++); 
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
        if (globalstatus && (!ix_update)) mssleep(DT);
    }
    PRINTF_DBG("Final termination of perform_requests :-)\n");std::cout<<std::flush;

} // end of function


#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
