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

void send_nonblocking(int owner, MPI_Request &r, double *ix, int TAG);
void recv_nonblocking(int owner, MPI_Request &r, double *result, int TAG);

void irespond_value(ReferenceContainer &REF, double ix, int owner, std::list<MPI_Request>::iterator &R);


//template<int DT, int TIMETOL, int BATCH>
//void answer_messages(ReferenceContainer &REF,
//                     int MYTHR,
//                     std::list<MPI_Request>::iterator &R,
//                     std::list<MPI_Request> &R_tot){};

template<int DT, int TIMETOL, int BATCH>
void answer_messages(ReferenceContainer &REF,
                     int MYTHR
//                     ,std::list<MPI_Request>::iterator &R,
//                     std::list<MPI_Request> &R_tot
                     ){

    std::list<MPI_Request> R_tot(10);
    std::list<MPI_Request>::iterator R = R_tot.begin();

// Potentially exetend the request strip :O

    // please add BATCH * (TIMETOL + 1) if N_remaining isnt enough
    int NDISPATCHED = 0;
    size_t d = std::distance(R, R_tot.end());
    int N_remaining = sizeof(d)/sizeof(MPI_Request);
    int N_SafeRequired = BATCH * (TIMETOL + 1) - N_remaining;
    if (N_SafeRequired > 0){
        for (int i=0; i<N_SafeRequired; i++){
            R_tot.push_back(MPI_Request());
        }
    }

    // lay the probes for all the procs
    int flagprobes[REF.p_ComHelper->WORLD_SIZE[MYTHR]] = {0};
    for (int i=0; i<REF.p_ComHelper->WORLD_SIZE[MYTHR]; i++){
        if (i != REF.p_ComHelper->WORLD_RANK[MYTHR]){
            MPI_Iprobe(i, 1, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
        }
    }


    // initialize auxiliary variables
    int ticks = 0;
    std::set<int> answered;
    answered.insert(REF.p_ComHelper->WORLD_RANK[MYTHR]);
    double ix;
    int i = 0;
    int status = 0;
    MPI_Request localreq;
    int status_localreq = 0;

    // start iterating
    while (ticks < TIMETOL){
        for (int j=0; j < BATCH; ++j) {
            if (flagprobes[i] == 1) { // a message is available! recieve it!
                printf("THIS IS AN ERROR: THERE SHOULD BE NO INCOMING MESSAGES!");
                // FLAG IS TRUE: TRY TO CAPTURE IT
                //
                // MPI_Recv(&ix, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // Blocking was not suited :-( So we do the following non-blocking :-)
                //
                MPI_Irecv(&ix,
                          1,//count
                          MPI_DOUBLE, // type
                          i, // destination
                          1, MPI_COMM_WORLD, &localreq);
                MPI_Request_get_status(localreq, &status_localreq, MPI_STATUS_IGNORE);
                if (status_localreq){ // IF INDEED THERE WAS  A MESSAGE, PROCEED.
                    if (R != R_tot.end()) {
                        irespond_value(REF, ix, i, R);
                        ++R;
                    }
                    flagprobes[i] = 0;
                    answered.insert(i);
                    NDISPATCHED++;
                    status_localreq = 0;
                } else { // CLEAN
                    flagprobes[i] = 0;
                    answered.insert(i); // it was answered, just not by us :-)
                }
            }
            if (i >= REF.p_ComHelper->WORLD_SIZE[MYTHR]) {
                i = 0;
            } else {
                ++i;
            }
            if (answered.size() == REF.p_ComHelper->WORLD_SIZE[MYTHR]){
                return; // we responded to all processes at least once, now we can leave :-)
                // this is a moment to note with clarity that if our graph has N nodes, P processes
                // and T threads devoted to answering operations (i.e. ~T/2 with T total threads),
                // then assuming that P-2 processes have 0 nodes, we have 1, and processor P-1 has
                // N-1 nodes that require our information, then we should have at least (N-1)
                // communications open with each processor. Here it is enforced in parallel
                // quite heuristically, but later in the singlestep_evolution code the all2all
                // is alternated with this answering role until all information is propagated.
            }
            if (answered.count(i) == 1) {
                    ++i;
            }
        }
        mssleep(DT);
        ++ticks;
//        std::cout << "Ticks are " << ticks << " out of " << TIMETOL << std::endl;
    }
    // leave
    //printf("Successfully dispatched %d requests\n", NDISPATCHED);
};

//template <int BATCH>
//void GetOneMsg(int ix,
//               ReferenceContainer &REF,
//               unsigned long N){};

template <int BATCH>
void GetOneMsg(int ix,
               ReferenceContainer &REF,
               unsigned long N){

//    std::cout << "Hi, welcome to 'GetOneMsg'... " << std::endl;
//    mssleep(500);

    if (BATCH<=0){error_report(min_batch_msg);} // guard: we need a batch size of at least 1 :-)


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
#pragma atomic
    auto itBeg = REF.p_ParHelper->data[ix].MissingA.begin();
#pragma atomic
    auto itEnd = REF.p_ParHelper->data[ix].MissingA.end();
    for (auto _it= itBeg; _it != itEnd; ++_it) {
        auto &thread = *_it;
        for (auto it =  thread.begin();it != thread.end(); it++){
            // retrieve data
            owner[QAvailable.front()] = std::get<1>(*it);
            their_vix[QAvailable.front()] = std::get<2>(*it);
            // send request for missing info
            send_nonblocking(owner[QAvailable.front()],
                                     requests_send[QAvailable.front()],
                                     & their_vix[QAvailable.front()], 1);
            recv_nonblocking(owner[QAvailable.front()],
                                     requests_recv[QAvailable.front()],
                                     & vval[QAvailable.front()], 3);
            results[QAvailable.front()] = std::make_tuple(0, // placeholder until we get the correct val
                                                          std::get<0>(*it),
                                                          owner[QAvailable.front()] * N + their_vix[QAvailable.front()]);
            QPend.push_back(QAvailable.front());
            QAvailable.pop();
            if (QAvailable.empty()) {
                waiting = true;
                while (waiting) {
                    auto i = QPend.begin();
                    while (i != QPend.end()) {
                        MPI_Request_get_status(requests_send[*i], &sstatus, MPI_STATUS_IGNORE);
                        MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                        if ((sstatus==1) && (rstatus==1)) {
                            std::get<0>(results[*i]) = vval[*i];
                            // No lock is needed: this index is entirely ours :-)
                            (*REF.p_IntHelper)[ix].ResultsPendProcess.push_back(results[*i]);
                            QAvailable.push(*i);
                            QPend.erase(i++);
                            waiting = false;
                        } else {
                            ++i;
                        }
                    }
                }
            }
        }
    }
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
}

template <int TIMETOL, int DT, int MAX_SUBTHR, int BATCH>
void GetAllMsgs(int NNodes,
                ReferenceContainer REF,
                unsigned long N,
                OpenMPHelper &O) {

//    std::cout << "Hi, welcome to 'GetAllMsgs'... " << std::endl;
//    mssleep(2000);

    // template the number of subthreads :-)
    if (MAX_SUBTHR<=0){error_report(min_subthread_msg);} // guard: we need at least 1 subth :-)
    //omp_set_dynamic(0);
    omp_set_num_threads(MAX_SUBTHR + 1);
    int N_SPAWNED = 0;
    int macro_threadn;
    macro_threadn = O.MY_THREAD_n;
    std::set<int> covered;
    long ix = - 1; // the index we will get, use and change from the global queue :-)!
//    std::cout << "Confirm we are at the face of the parallel region..." << std::endl;
//    mssleep(2000);

#pragma omp parallel firstprivate(REF, macro_threadn, NNodes, N)
{
    if (omp_get_num_threads()==1) {
        //PRINTF_DBG("One subthread is not enough; recursive restart.");
//        std::cout << "There is only one thread, so I will answer some messages and then call again itself hoping that more threads are spawned :o"<< std::endl;
        //answer_messages<DT, TIMETOL, BATCH>(REF, macro_threadn, p_REQUESTLIST, REQUESTLIST);
    }
    //} else if ((omp_get_thread_num() == 0) && (omp_get_num_threads()>1)) {
    if ((omp_get_thread_num() == 0) && (omp_get_num_threads()>1)) {
        bool spawned_max_capacity = true;
        bool isempty = false;
        bool globalstatus = true;
        int has_been_consumed=0;
        int atomic_helper;
#pragma omp atomic read // check if this loop still makes sense
        atomic_helper = *(REF.p_TOT);
        globalstatus = (atomic_helper < NNodes);
        while (globalstatus) {
#pragma omp critical
{
            isempty = REF.p_CHECKED->empty();
}
#pragma omp read
            atomic_helper = N_SPAWNED;
            spawned_max_capacity = (atomic_helper == MAX_SUBTHR);
            if ((spawned_max_capacity) || (isempty)) {
                answer_messages<DT, TIMETOL, BATCH>(REF, macro_threadn);
            } else if (!isempty) { // update our index so
                // in the elapsed time, it could happen that the only available
                // thread has just become busy and now you would be overwriting a value
                // i.e. (1) a new "ix" and an idle thread exist
                //      (2) you interpret the idle thread as the new "ix" having been consumed,
                //          and as the global queue isn't empty you will overwrite ix with that value
                //      (3) you could overwrite it before a thread is able to consume it, so instead
                //          of just using a queue for this we patch this problem by double-checking
                //          with a (locally) shared set.
                if (VERBOSE) {
                    // use PRINTF_DBG()
                    std::cout << " MASTER: capacity is not max, and an ix appears available... " << std::endl;
                }
#pragma omp critical
{
                    if ((ix==-1) || (covered.count(ix) == 1)) {
                        ix = REF.p_CHECKED->front();
                        REF.p_CHECKED->pop();
//                        std::cout << " MASTER: index update was effective" << std::endl;
                    }
}
                }
#pragma omp atomic read // check if this loop still makes sense
                atomic_helper = *(REF.p_TOT);
                globalstatus = (atomic_helper < NNodes);
        }
        // Please spend some prudential time answering messages before exiting to join the communal integration :O
        //PRINTF_DBG("Skipped answering messages :-)\n");
        //answer_messages<DT, TIMETOL, BATCH>(REF, macro_threadn, p_REQUESTLIST, REQUESTLIST);
    } else if (omp_get_num_threads()>1) {
        bool globalstatus = true;
        int local_ix;
        int old_ix = -1;
#pragma omp atomic read
        local_ix = ix;
        bool ix_update = false;
        int atomic_helper;
        //if (omp_get_thread_num() == 1){ix_update = true;} // only thread=1 grabs the initial index
#pragma omp atomic read
        atomic_helper = *(REF.p_TOT);
        globalstatus = (atomic_helper < NNodes); // just in case: update global status. At first it should
        //   not be necessary but it could occur in relatively small systems and in single-processor cases
//        std::cout << " global_status and ix_update are initially: " << globalstatus << " and " << ix_update << std::endl;
        while (globalstatus) {
            //mssleep(5000);
            //std::cout << " subworker " << omp_get_thread_num() << " of worker " << macro_threadn << " has completed one lap" << std::endl;
            if (ix_update) { // Perform the main task: ask for the required information!
#pragma omp atomic update
                ++ N_SPAWNED;
//                std::cout << "The other thread will go try get one msg... " << std::endl;
//                mssleep(1000);
                //std::cout << "[DEBUG] about to call GetOneMsg" << std::endl;
                //GetOneMsg<BATCH>(local_ix, REF, N);
                //PRINTF_DBG("Skipped sending messages :-)\n");
#pragma omp atomic update
                ++ *(REF.p_TOT);
                // Also communicate to other threads that there is one new item pending integration :o
#pragma critical
{
                REF.p_READY_FOR_INTEGRATION->push(local_ix);

}
#pragma omp atomic update
                -- N_SPAWNED; // let the thread=0 know you may become idle.
                ix_update = false;
//                std::cout << " subworker finished processing data via GetOneMsg " << ix << std::endl;
            }
#pragma omp atomic read // recompute the index: maybe a new one is available!
            local_ix = ix;
            if (local_ix != old_ix) { // if there is a new one, please make sure it hasnt been
#pragma omp critical                // grabbed by some other thread: attempt to claim sole ownership
{                                   // in a critically-protected block ;-).
                if (!covered.count(ix)) {
                    ix_update = true;
                    covered.insert(ix); // effectively claim ownership
//                    std::cout << " subworker claimed ownership of ix " << ix << std::endl;
                }
}               // outside the critical block:
                if (!ix_update){old_ix=local_ix;} // let it go: it has already been claimed
            }
#pragma omp atomic read // check if this loop still makes sense
            atomic_helper = *(REF.p_TOT);
            globalstatus = (atomic_helper < NNodes);
        }
        // Please spend some prudential time answering messages before joining the communal integration :O
        //answer_messages<DT, TIMETOL, BATCH>(REF, macro_threadn, p_REQUESTLIST, REQUESTLIST);
        //PRINTF_DBG("Skipped answering messages :-)\n");
    } // this was all for the threads!= 0
} // end of parallel region
} // end of function





void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
