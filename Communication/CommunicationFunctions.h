//
// Created by m4zz31 on 5/11/21.
//

#ifndef CPPPROJCT_COMMUNICATIONFUNCTIONS_H
#define CPPPROJCT_COMMUNICATIONFUNCTIONS_H

#include <boost/serialization/string.hpp>
#include "../GraphClasses/GeneralGraph.h"
#include "../GraphClasses/GraphFunctions.h"
#include "../Utils/HelperClasses.h"
#include "../Utils/error.h"
#include "../Utils/global_standard_messages.h"
#include <set>

void send_nonblocking_request(int &owner, MPI_Request &r, double *ix, int TAG);
void recv_nonblocking_request(int &owner, MPI_Request &r, double *result, int TAG);

template <int BATCH>
void GetOneMsg(int ix,
               ReferenceContainer &REF,
               unsigned long N){};

//template <int BATCH>
//void GetOneMsg_BIS(int ix,
//               ReferenceContainer &REF,
//               unsigned long N){
//
//    if (BATCH<=0){error_report(min_batch_msg);} // guard: we need a batch size of at least 1 :-)
//
//
//    bool waiting = false;
//    int rstatus = 0, sstatus = 0;
//    MPI_Request requests_send[BATCH], requests_recv[BATCH];
//    InfoVecElem results[BATCH];
//    std::queue<int> QAvailable;
//    std::list<int> QPend;
//    for (int i = 0; i < BATCH; ++i){
//        QAvailable.push(i);
//    }
//
//    double vval[BATCH], their_vix[BATCH];
//    int owner[BATCH];
//
//    // Collect missing nodes
//    for (auto & thread : P.data[ix].MissingA) {
//        for (auto it =  thread.begin();it != thread.end(); it++){
//            // retrieve data
//            owner[QAvailable.front()] = std::get<1>(*it);
//            their_vix[QAvailable.front()] = std::get<2>(*it);
//            // send request for missing info
//            send_nonblocking_request(owner[QAvailable.front()],
//                                     requests_send[QAvailable.front()],
//                                     & their_vix[QAvailable.front()], 1);
//            recv_nonblocking_request(owner[QAvailable.front()],
//                                     requests_recv[QAvailable.front()],
//                                     & vval[QAvailable.front()], 3);
//            results[QAvailable.front()] = std::make_tuple(0, // placeholder until we get the correct val
//                                                          std::get<0>(*it),
//                                                          owner[QAvailable.front()] * N + their_vix[QAvailable.front()]);
//            QPend.push_back(QAvailable.front());
//            QAvailable.pop();
//            if (QAvailable.empty()) {
//                waiting = true;
//                while (waiting) {
//                    auto i = QPend.begin();
//                    while (i != QPend.end()) {
//                        MPI_Request_get_status(requests_send[*i], &sstatus, MPI_STATUS_IGNORE);
//                        MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
//                        if ((sstatus==1) && (rstatus==1)) {
//                            std::get<0>(results[*i]) = vval[*i];
//#pragma omp atomic
//                            I[ix].ResultsPendProcess.push_back(results[*i]);
//                            QAvailable.push(*i);
//                            QPend.erase(i++);
//                            waiting = false;
//                        } else {
//                            ++i;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    // Collect missing nodes+edges
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
//}

template <int TIMETOL, int DT, int MAX_SUBTHR, int BATCH>
void GetAllMsgs(int NNodes,
                ReferenceContainer REF,
                unsigned long N,
                OpenMPHelper &O) {

    // template the number of subthreads :-)
    if (MAX_SUBTHR<=0){error_report(min_subthread_msg);} // guard: we need at least 1 subth :-)
    omp_set_num_threads(MAX_SUBTHR + 1);
    int N_SPAWNED = 0;
    int macro_threadn;
    macro_threadn = O.MY_THREAD_n;
    std::set<int> covered;
    long ix; // the index we will get, use and change from the global queue :-)!

#pragma omp parallel firstprivate(REF, macro_threadn)
{
    if (omp_get_thread_num() == 0) {
        bool spawned_max_capacity = true;
        bool isempty = false;
        bool globalstatus = true;
        int has_been_consumed=0;
        int atomic_helper;
#pragma omp atomic read // check if this loop still makes sense
        atomic_helper = *(REF.p_TOT);
        globalstatus = (atomic_helper < NNodes);
        while (globalstatus) {
            mssleep(500);
            std::cout << " subworker " << omp_get_thread_num() << " of worker " << macro_threadn << " has completed one lap" << std::endl;
#pragma omp critical
{
            isempty = REF.p_CHECKED->empty();
}
#pragma omp read
            atomic_helper = N_SPAWNED;
            spawned_max_capacity = (atomic_helper == MAX_SUBTHR);
            if ((spawned_max_capacity) || (isempty)) {
                // spend some time answering messages :-0
                mssleep(DT * TIMETOL);
            } else if (!isempty) { // update our index so
                // in the elapsed time, it could happen that the only available
                // thread has just become busy and now you would be overwriting a value
                // i.e. (1) a new "ix" and an idle thread exist
                //      (2) you interpret the idle thread as the new "ix" having been consumed,
                //          and as the global queue isn't empty you will overwrite ix with that value
                //      (3) you could overwrite it before a thread is able to consume it, so instead
                //          of just using a queue for this we patch this problem by double-checking
                //          with a (locally) shared set.
#pragma omp critical
{
                    if (covered.count(ix) == 1) {
                        ix = REF.p_CHECKED->front();
                        REF.p_CHECKED->pop();
                    }
}
                }
#pragma omp atomic read // check if this loop still makes sense
                atomic_helper = *(REF.p_TOT);
                globalstatus = (atomic_helper < NNodes);
        }
    }
    else {
        bool globalstatus = true;
        int local_ix;
#pragma omp atomic read
        local_ix = ix;
        int old_ix = local_ix;
        bool ix_update = false;
        int atomic_helper;
        if (omp_get_thread_num() == 1){ix_update = true;} // only thread=1 grabs the initial index
#pragma omp atomic read
        atomic_helper = *(REF.p_TOT);
        globalstatus = (atomic_helper < NNodes); // just in case: update global status. At first it should
        //   not be necessary but it could occur in relatively small systems and in single-processor cases
        while (globalstatus) {
            mssleep(500);
            std::cout << " subworker " << omp_get_thread_num() << " of worker " << macro_threadn << " has completed one lap" << std::endl;
            if (ix_update) { // Perform the main task: ask for the required information!
#pragma omp atomic update
                ++ N_SPAWNED;
                GetOneMsg<BATCH>(local_ix, REF, N);
#pragma omp atomic update
                ++ *(REF.p_TOT);
#pragma omp atomic update
                -- N_SPAWNED; // let the thread=0 know you may become idle.
                ix_update = false;
            }
            int local_ix;
#pragma omp atomic read // recompute the index: maybe a new one is available!
            local_ix = ix;
            if (local_ix != old_ix) { // if there is a new one, please make sure it hasnt been
#pragma omp critical                // grabbed by some other thread: attempt to claim sole ownership
{                                   // in a critically-protected block ;-).
                if (!covered.count(ix)) {
                    ix_update = true;
                    covered.insert(ix); // effectively claim ownership
                }
}               // outside the critical block:
                if (!ix_update){old_ix=local_ix;} // let it go: it has already been claimed
            }
#pragma omp atomic read // check if this loop still makes sense
            atomic_helper = *(REF.p_TOT);
            globalstatus = (atomic_helper < NNodes);
        }
    } // this was all for the threads!= 0
} // end of parallel region

    // Please spend some prudential time answering messages before joining the communal integration :O
    mssleep(TIMETOL * DT);

} // end of function





void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
