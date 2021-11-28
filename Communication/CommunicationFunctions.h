//
// Created by m4zz31 on 5/11/21.
//

#ifndef CPPPROJCT_COMMUNICATIONFUNCTIONS_H
#define CPPPROJCT_COMMUNICATIONFUNCTIONS_H

#include <set>
#include <iterator>
#include <random>
#include <type_traits>

#include "../macros/macros.h"

#include "../GraphClasses/GeneralGraph.h"
#include "../GraphClasses/GraphFunctions.h"

#include "../Utils/HelperClasses.h"
#include "../Utils/error.h"
#include "../Utils/global_standard_messages.h"
#include "../Utils/msleep.h"



void build_answer_nodes(double &answer, ReferenceContainer &REF, double ix, int owner, int &MyNProc);
void build_answer_edges(double * answer, ReferenceContainer &REF, double * ix, int owner, int &MyNProc);


template <typename SpecificRequestObject>
class RequestObject : public SpecificRequestObject {
    public:
        RequestObject(int type): SpecificRequestObject(type){};
        RequestObject(){};
};


template <int field>
class FieldRequestObject {
public:
    int recvLength = -1;
    int sendLength = -1;
    int recvTag = -1;
    int sendTag = -1;
    bool recvInt = true;
    bool sendDouble = true;
    // Reception type
    typedef typename std::conditional< (field >= 0), int, double>::type buffer_type;
    typedef typename std::conditional< (field >= -2), double, double>::type answer_type;
    FieldRequestObject();
    // Diffrent SendTag builders
    void buildSendTag(int * data);
    void buildSendTag(double * data);
    // Different RecvTag builders
    void buildRecvTag(int * data){};
    void buildRecvTag(double * data){};
    void buildRecvTag();
    // integer buffer
    void computeAnswer(ReferenceContainer &REF, int * buffer, double * answer, MPI_Status &S, int MYPROC);
    void computeReady(ReferenceContainer &REF, int * buffer, bool &isReady);
    // double buffer
    void computeReady(ReferenceContainer &REF, double * buffer, bool &isReady);
    void computeAnswer(ReferenceContainer &REF, double * buffer, double * answer, MPI_Status &S, int MYPROC);
};

template <int field>
void FieldRequestObject<field>::buildSendTag(int * data){
    sendTag = *data;
};

template <int field>
void FieldRequestObject<field>::buildSendTag(double * data){
    if (field == -1){
        sendTag = (int) (*data);
    } else if (field == -2) {
        sendTag = (int) ((int) (*(data + 1)) + (int) OFFSET);
    }
};


template <int field>
void FieldRequestObject<field>::buildRecvTag(){
    if (field == 0){
        recvTag = K1_REQUEST;
        return;
    } else if (field == 1) {
        recvTag = K2_REQUEST;
        return;
    } else if (field == 2) {
        recvTag = K3_REQUEST;
        return;
    } else if (field == -1) {
        recvTag = VERTEXVAL_REQUEST_FLAG;
        return;
    } else if (field == -2){
        recvTag = EDGEVAL_REQUEST_FLAG;
        return;
    }
};



template <int field>
FieldRequestObject<field>::FieldRequestObject(){
    if (field == 0){
        // Answering requests for the field term 1
        recvLength = 1;
        sendLength = 1;

    } else if (field == 1) {
        // Answering requests for the field term 2
        recvLength = 1;
        sendLength = 1;

    } else if (field == 2) {
        // Answering requests for the field term 3
        recvLength = 1;
        sendLength = 1;

    } else if (field == -1) {
        // Answering requests for the node values
        recvInt = false;
        sendDouble = true;
        recvLength = 1;
        sendLength = 1;

    } else if (field == -2) {
        // Answering requests for the edge values
        recvLength = 2;
        sendLength = 2;
        recvInt = false;
        sendDouble = true;
    }
}

template <int field>
void FieldRequestObject<field>::computeAnswer(ReferenceContainer &REF,
                                       int * buffer,
                                       double * answer,
                                       MPI_Status &S,
                                       int MYPROC){
    if (field == 0){
#pragma omp atomic read
        (*answer) = REF.p_LayHelper->data[*(buffer)].RK1[0];
    } else if (field == 1) {
#pragma omp atomic read
        (*answer) = REF.p_LayHelper->data[*(buffer)].RK2[0];
    } else if (field == 2) {
#pragma omp atomic read
        (*answer) = REF.p_LayHelper->data[*(buffer)].RK3[0];
    }
};


template <int field>
void FieldRequestObject<field>::computeAnswer(ReferenceContainer &REF,
                                       double * buffer,
                                       double * answer,
                                       MPI_Status &S,
                                       int MYPROC){
    if (field == -1) {
        build_answer_nodes(*answer, REF, *buffer, S.MPI_SOURCE, MYPROC);
    } else if (field == -2) {
        build_answer_edges(answer, REF, buffer, S.MPI_SOURCE, MYPROC);
    }
};



template <int field>
void FieldRequestObject<field>::computeReady(ReferenceContainer &REF, double * buffer, bool &isReady){
    if (field == -1) {
        isReady = true;
        return;
    } else if (field == -2) {
        isReady = true;
        return;
    }
};

template <int field>
void FieldRequestObject<field>::computeReady(ReferenceContainer &REF, int * buffer, bool &isReady){
    if (field == 0){
#pragma omp atomic read
        isReady = REF.p_LayHelper->data[*(buffer)].RK1_status;
        return;
    } else if (field == 1) {
#pragma omp atomic read
        isReady = REF.p_LayHelper->data[*(buffer)].RK2_status;
        return;
    } else if (field == 2) {
#pragma omp atomic read
        isReady = REF.p_LayHelper->data[*(buffer)].RK3_status;
        return;
    }
};




template<int DT, int TIMETOL, int BATCH, typename RequestClass>
void generic_answer_requests(ReferenceContainer &REF, int MYTHR, RequestClass ReqObj){
    int flag = 0;
    int TRIES=0;
    int MYPROC = REF.p_ComHelper->WORLD_RANK[MYTHR];
    MPI_Request R;
    MPI_Message M;
    MPI_Status S;
    int t=0;
    bool firstlap = true;

    // Should we receive an int or a double
    typename decltype(ReqObj)::buffer_type buffer[ReqObj.recvLength];
    typename decltype(ReqObj)::answer_type answer[ReqObj.sendLength];

    ReqObj.buildRecvTag();
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
        PRINTF_DBG("About to probe! recvTag is %d\n", ReqObj.recvTag);
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
            PRINTF_DBG("[GAR] About to receive! recvInt is %d, and recvLength is %d\n", ReqObj.recvInt, ReqObj.recvLength);
            if (ReqObj.recvInt){
                MPI_Mrecv(&buffer,
                          ReqObj.recvLength,
                          MPI_INT,
                          &M, &S);
            } else {
                MPI_Mrecv(&buffer,
                          ReqObj.recvLength,
                          MPI_DOUBLE,
                          &M, &S);
            }
            PRINTF_DBG("[GAR] Received successfully!\n");
            bool isReady = false;
#pragma omp critical
{
            ReqObj.computeReady(REF, &buffer[0], isReady);
};
            while (!isReady) {
#pragma omp critical
{
                ReqObj.computeReady(REF, &buffer[0], isReady);
};
                PRINTF_DBG("[GAR] not ready yet!\n");
                mssleep(DT);
            }
            ReqObj.computeAnswer(REF, &buffer[0], &answer[0], S, MYPROC);
            ReqObj.buildSendTag(&buffer[0]);
            PRINTF_DBG("[GAR] About to send! sendDouble is %d, ReqObj.sendLength is %d, and ReqObj.sendTag is %d\n", ReqObj.sendDouble, ReqObj.sendLength, ReqObj.sendTag);
            PRINTF_DBG("[GAR] The sent message('s first element) will be %f\n", answer[0]);
            if (ReqObj.sendDouble){
                MPI_Ssend(&answer[0],
                          ReqObj.sendLength,
                          MPI_DOUBLE,
                          S.MPI_SOURCE,
                          ReqObj.sendTag,
                          MPI_COMM_WORLD);
            } else {
                MPI_Ssend(&answer[0],
                          ReqObj.sendLength,
                          MPI_INT,
                          S.MPI_SOURCE,
                          ReqObj.sendTag,
                          MPI_COMM_WORLD);
            }
        }
    }
}















template<int DT, int TIMETOL, int BATCH>
void perform_field_requests(ReferenceContainer &REF,int MYPROC, int fieldOrder,std::queue<long> * queue){

    int ASKING_TAGS[4] = {VERTEXVAL_REQUEST_FLAG, K1_REQUEST, K2_REQUEST, K3_REQUEST};
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
        } else if (fieldOrder==1) {
            assert((*REF.p_IntHelper)[ix].ixMap.size() == (REF.p_LayHelper->data[ix].RK2.size()-1));
            L = (*REF.p_IntHelper)[ix].ixMap.size() + 1;
        } else {
            printf("[FATAL] field order requested to perform_field_requests does not exist!\n");std::cout<<std::flush;
            exit(1);
        }

        for (int i=0; i<L-1 ; ++i){
            double recvBuffer;
            if (((int) std::get<2>((*REF.p_IntHelper)[ix].ixMap[i])) == MYPROC) {
                unsigned long owner = std::get<1>((*REF.p_IntHelper)[ix].ixMap[i]);
                bool isReadyYet = false;
                if (fieldOrder == 2) {
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
#pragma omp atomic read
                    isReadyYet = REF.p_LayHelper->data[owner].RK3_status;
                    while (!isReadyYet) {
                        mssleep(DT);
#pragma omp atomic read
                        isReadyYet = REF.p_LayHelper->data[owner].RK3_status;
                    }
#pragma omp atomic read
                    recvBuffer = REF.p_LayHelper->data[owner].RK3[0];
                } else if (fieldOrder == 1) {
#pragma omp atomic read
                    recvBuffer = (*REF.p_IntHelper)[owner].centralValue;
                };
            } else {

                // This is kinda dumb, but the 'else' clause is for compatibility.
                if (fieldOrder != 1){
                    int sendBuffer = (int) std::get<1>((*REF.p_IntHelper)[ix].ixMap[i]);
                    PRINTF_DBG("[pfr]!=1 About to send! sendbuffer says %d... asking tag is %d\n", sendBuffer, ASKING_TAGS[fieldOrder-1]);
                    MPI_Ssend(&sendBuffer,
                              1,
                              MPI_INT,
                              (int) std::get<2>((*REF.p_IntHelper)[ix].ixMap[i]),
                              ASKING_TAGS[fieldOrder-1],
                              MPI_COMM_WORLD);
                    PRINTF_DBG("Sent successfully!\n");
                } else {
                    double sendBuffer = (double) std::get<1>((*REF.p_IntHelper)[ix].ixMap[i]);
                    PRINTF_DBG("[pfr]==1 About to send! sendbuffer says %f... asking tag is %d\n", sendBuffer, ASKING_TAGS[fieldOrder-1]);
                    MPI_Ssend(&sendBuffer,
                              1,
                              MPI_DOUBLE,
                              (int) std::get<2>((*REF.p_IntHelper)[ix].ixMap[i]),
                              ASKING_TAGS[fieldOrder-1],
                              MPI_COMM_WORLD);
                    PRINTF_DBG("Sent successfully!\n");
                }
                PRINTF_DBG("[pfr] About to receive! tag is %d\n", (int) std::get<1>((*REF.p_IntHelper)[ix].ixMap[i]));
                MPI_Recv(&recvBuffer,
                         1,
                         MPI_DOUBLE,
                         (int) std::get<2>((*REF.p_IntHelper)[ix].ixMap[i]),
                         (int) std::get<1>((*REF.p_IntHelper)[ix].ixMap[i]),
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
            } else if (fieldOrder == 1) {
                (*REF.p_IntHelper)[ix].neighborValues[i] =  recvBuffer;
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
