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

void irespond_value(ReferenceContainer &REF, double ix, int owner, MPI_Request & R, int MyNProc);

void recv_blocking(int owner, MPI_Request &r, double &result, int TAG);


template<int DT, int TIMETOL, int BATCH>
void answer_messages(ReferenceContainer &REF,int MYTHR) {

    int MYPROC = REF.p_ComHelper->WORLD_RANK[MYTHR];
    int NPROCS = REF.p_ComHelper->WORLD_SIZE[MYTHR];
    std::vector<MPI_Request> R_tot;
    std::vector<MPI_Request> R_send;
//    R_tot.push_back(MPI_Request());
//    R_send.push_back(MPI_Request());
    int NDISPATCHED = 0;

    // lay the probes for all the p rocs
    int flagprobes[NPROCS];
    int probe_status = 1;
    int statusreq_status = 1;
    for (int i = 0; i < NPROCS; i++) {
        flagprobes[i] = 0;
        if (i != MYPROC) {
            probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
            while (probe_status!=0){
                probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
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
    double ix;
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

            // re-probe this index ;-)
//            probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//            while (probe_status != 0){
//                probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//            }
            // fake-probe this index :-<
            // just to debug: in the future we should 'fake probe' just once in a while ;-)
            //flagprobes[i] = 1;

            //printf("running inner loop\n");
            if (flagprobes[i] == 1) { // a message is available! recieve it!
                PRINTF_DBG("Entered flagprobes\n");
                R_tot.push_back(MPI_Request());

                if (true) {
                    recv_nonblocking(i, R_tot[R_tot.size() - 1], ix, VERTEXVAL_REQUEST_FLAG);

                    printf("|%d %d %d|\n",flagprobes[0] , flagprobes[1] , flagprobes[2]);std::cout<<std::flush;

// NONINVASIVE APPROACH (A)
//                statusreq_status = MPI_Request_get_status(R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
//                while (statusreq_status != 0){
//                    mssleep(10); // HYPERPARAM
//                    printf("I am not being able to get the status!\n");
//                    statusreq_status = MPI_Request_get_status(R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
//                }
// INVASIVE APPROACH: (B)
                    statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
                    while (statusreq_status != 0) {
                        mssleep(10); // HYPERPARAM
                        statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
                    }
                } else {
                    recv_blocking(i, R_tot[R_tot.size() - 1], ix, VERTEXVAL_REQUEST_FLAG);
                    status_localreq = 1;
                }

                if (status_localreq == 1) { // If we were first to capture the message, proceed.
                    R_send.push_back(MPI_Request());
                    //printf("We effectively captured a vertex info request :-)\n");
                    irespond_value(REF, ix, i, R_send[R_send.size() - 1], MYPROC);
                    ++NTOT;
                    //printf("We effectively answered asynchronously a vertex info request :-)\n");
                    NDISPATCHED++;
                } else {
                    //printf(
                    //        "We were faced with a probe that indicated an incoming message but we couldnt capture it :o\n");
                    //std::cout << std::flush;
                }

                // ONLY NECESSARY IF WE ARE IGNORING THE PROBE AND JUST TRYING JUST IN CASE
                statusFree = MPI_Request_free(&R_tot[R_tot.size() - 1]);
                //if (statusFree != 0)  printf("answer_messages failed to free one req :O it was%d\n", statusFree);
                status_localreq = 0;
                flagprobes[i] = 0;

                // re-probe it
//                probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//                while (probe_status != 0){
//                    probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//                }

                // --------------------USELESS... CROP IT..............................
//                answered.insert(i); // it was answered, just not by us :-)
//            }
//            if (answered.size() == NPROCS) {
//                // we should just refresh answered :-)
//                answered = std::set<int>();
//                answered.insert(MYPROC);
//                // and reseed the probe :-)
//                for (int i = 0; i < NPROCS; i++) {
//                    flagprobes[i] = 0;
//                    if (i != MYPROC) {
//                        MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
//                    }
//                }
//                // or maybe just return now, i.e. earlier.
//                i = -1;

            // USELES-------REMOVE-----
            }

            //PRINTF_DBG("Arrived to C\n");
            ++i;
            if (i == MYPROC) ++i;
            if (i >= NPROCS) {
                //PRINTF_DBG("Entered D\n");
                if (MYPROC!=0) {
                    i = 0;
                } else {
                    i = 1;
                }
            }

            // USELES-------REMOVE-----
//            while (answered.count(i) == 1) {
//                //printf("Entered E\n");
//                ++i;
//                if (i == MYPROC) ++i;
//                if (i >= NPROCS) {
//                    i = 0;
//                }
//            }
            // USELES-------REMOVE-----


        } // END OF THE BATCH... (blocking) wait for the unmatched requests :-)
        if ((R_tot.size() > 0) || (R_send.size() > 0)) {
            PRINTF_DBG("ENDing answer_messages. Caught (R_tot)=%d requests, Answered (R_send)=%d.\n", R_tot.size(),
                       R_send.size());
            std::cout << std::flush;
            int status_of_getstatus = 1;
            for (int i = 0; i < R_send.size(); ++i) { // we wait for our answers to arrive before descoping

                // GET STATUS VERSIE
//                status_of_getstatus = MPI_Request_get_status(R_send[i], &status_localreq, MPI_STATUS_IGNORE);
//                while (status_of_getstatus != 0) {
//                    status_of_getstatus = MPI_Request_get_status(R_send[i], &status_localreq, MPI_STATUS_IGNORE);
//                    mssleep(5); // HYPERPARAM
//                }
                // TEST VERSIE
                status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
                while (status_of_getstatus != 0) {
                    printf("Failing to test a request, its returning %d\n", status_of_getstatus);
                    status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
                    mssleep(5); // HYPERPARAM
                }

                if (status_localreq == 0) {
                    for (int jk=0;jk<100;++jk) printf("We are currently waiting for a send to be completed :-(\n");
                    std::cout << std::flush;
                    MPI_Wait(&R_send[i], MPI_STATUS_IGNORE);
                }
                //destroyRequest(R_send[i], NERR);
                PRINTF_DBG("Successfully waited the request to be completed\n");//, now we will attempt to free it.\n");
                // ATTEMPTING TO FREE IS POINTLESS NOW, AS IT WAS A BLOCKING WAIT :-)
//                std::cout << std::flush;
//                statusFree = MPI_Request_free(&R_send[i]);
//                PRINTF_DBG("Successfully finished answering that one request! ^.^ status of free: %d\n", statusFree);
//                std::cout << std::flush;
            }
        }
        std::vector<MPI_Request> R_tot;
        std::vector<MPI_Request> R_send;
        mssleep(DT);
        ++ticks;
    }
    //printf("ENDed answering\n");
};







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
    long current_ix;

    bool all_sent_locally = false;
    int tot_locals = 0;
    int sent_locally = 0;

    // insert here required initialization
    std::list<long> all_indexes_seen;
    if (BATCH<=0)error_report(min_batch_msg); // guard: we need a batch size of at least 1 :-)
    bool waiting = false;
    bool was_last_appearance = false;
    int special_index = -1;
    int rstatus = 0, sstatus = 0;
    MPI_Request requests_send[BATCH], requests_recv[BATCH];
    InfoVecElem results[BATCH];
    std::queue<int> QAvailable;
    std::list<int> QPend;
    for (int i = 0; i < BATCH; ++i){
        QAvailable.push(i);
    }
    double vval[BATCH], their_vix[BATCH];
    int ixlist[BATCH]; // this is our new helper in order to process various indexes asynchronously (=^-|)-|-/
    int owner[BATCH];


#pragma omp atomic read
    atomic_helper = *(REF.p_TOT);
    globalstatus = (atomic_helper < NNodes);

    while (globalstatus) {
        // Perform the main task: ask for the required information!
        if (ix_update) {
//            printf("About to ask for the process of local elem ix: %d\n",ix);
            auto itBeg = REF.p_ParHelper->data[ix].MissingA.begin();
            auto itEnd = REF.p_ParHelper->data[ix].MissingA.end();
            tot_locals  = 0;
            sent_locally = 0;
            total_processed = 0;

            // race-condition unfriendly?
            for (auto _it= itBeg; _it != itEnd; ++_it) {
                auto &thread = *_it;
                tot_locals += itEnd - itBeg;

//                std::cout << "A" << std::endl;
//                std::cout << std::flush; // DEBUGGING
                // DEBUG SECTION
                for (auto it =  thread.begin();it != thread.end(); it++){

                    // retrieve data
                    owner[QAvailable.front()] = std::get<1>(*it);
                    their_vix[QAvailable.front()] = std::get<2>(*it);

//                    std::cout << "B" << std::endl;
//                    std::cout << std::flush; // DEBUGGING
                    // Debug! :-)
                    printf("Sending and recieving (nonblocking): asking %d for vertex w index: %f \n",
                           owner[QAvailable.front()], their_vix[QAvailable.front()]);
                    std::cout  << std::flush; // DEBUGGING

                    // Send request for missing info
                    send_nonblocking(owner[QAvailable.front()], // send a request (tag=-1) for some ix's vertex value
                                     requests_send[QAvailable.front()],
                                     their_vix[QAvailable.front()], VERTEXVAL_REQUEST_FLAG);

                    // Recieve the answer for the requested info, using tag = index
                    recv_nonblocking(owner[QAvailable.front()], // recieve the respone (tag=ix) for that  ix
                                     requests_recv[QAvailable.front()], // flag = (int) index :-)
                                     vval[QAvailable.front()], (int) their_vix[QAvailable.front()]);

                    // Store the result in our BATCH-sized temporal container:
                    // we are only lacking the exact vertex value which should be returned from the owner
                    results[QAvailable.front()] = std::make_tuple(0, // placeholder until we get the correct val
                                                                  std::get<0>(*it),
                                                                  owner[QAvailable.front()] * N + their_vix[QAvailable.front()]);

                    PRINTF_DBG("A\n");
                    std::cout << std::flush;
                    // Store the temporal results index that is requiring an answer
                    QPend.push_back(QAvailable.front());
                    PRINTF_DBG("B\n");
                    std::cout << std::flush;
                    // add to ixlist[QAvailable.front()] the ix
                    // so that if we dont get the answer in this ix then we can do it in the future
                    ixlist[QAvailable.front()] = ix;
                    PRINTF_DBG("C\n");
                    std::cout << std::flush;
                    // Remove that last used element from QAvailable ;-)
                    QAvailable.pop();
                    PRINTF_DBG("queue size is %d\n", QAvailable.size());
                    std::cout << std::flush;
                    // If QAvailable is empty, devote ourselves to asynchronously waiting for
                    // answers to arrive. This indirectly means that our QPend has reached length BATCH.
                    if (QAvailable.empty()) {
                        waiting = true;
                        PRINTF_DBG("it IS FINALLY EMPTY\n");
                        std::cout << std::flush;

                        while (waiting) {
                            sstatus = 0;
                            rstatus = 0;
                            printf("We are waiting to be responded :-( Qpends size is %d\n", QPend.size());
                            std::cout << std::flush;
                            //mssleep(150);
                            // Iterate through the pending indexes to see if one has been answered
                            auto i = QPend.begin();
                            auto status_sstatus=1, status_rstatus=1;
                            while (i != QPend.end()) {
                                // Check if the request of index *i has been recieved and if the answer
                                // has or not already arrived.
                                status_sstatus = MPI_Request_get_status(requests_send[*i], &sstatus, MPI_STATUS_IGNORE);
                                while (status_sstatus!=0){
                                    printf("fetching the status_sstatus failed :0 it yielded %d\n", status_sstatus);std::cout<<std::flush;
                                    status_sstatus = MPI_Request_get_status(requests_send[*i], &sstatus, MPI_STATUS_IGNORE);
                                }
                                status_rstatus = MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                                while (status_rstatus!=0){
                                    printf("fetching the status_rstatus failed it yielded %d:0\n", status_rstatus);std::cout<<std::flush;
                                    status_rstatus = MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                                }
                                // If it both arrived to them and was answered and came back, then
                                // store that value and free one space in QPend while inserting the
                                // new freed index into QAvailable.
                                printf("sstatus = %d, rstatus = %d, status_sstatus =%d, status_rstatus=%d\n", sstatus, rstatus,status_sstatus,status_rstatus);
                                std::cout << std::flush;
                                if ((sstatus==1) && (rstatus==1)) {
                                    // We keep track of the locally successfully sent and recieved
                                    // requests as to at last check if, for this index, we were able
                                    // to settle the issue and its ready for integration or not :-)
                                    std::get<0>(results[*i]) = vval[*i];
                                    special_index = ixlist[*i];
                                    (*REF.p_IntHelper)[special_index].ResultsPendProcess.push_back(results[*i]);
                                    PRINTF_DBG("Recieved a response to our request: val %f for ix: %d\n", vval[*i], special_index);
                                    std::cout << std::flush;
                                    // Set the index "*i" of the batch-sized containers as available
                                    QAvailable.push(*i);
                                    // erase the waiting clause
                                    waiting = false;
                                    // Regresh the MPI Request objects ;-)
                                    MPI_Request_free(&requests_send[*i]);
                                    MPI_Request_free(&requests_recv[*i]);
                                    requests_send[*i] = MPI_Request();
                                    requests_recv[*i] = MPI_Request();
                                    // QPend is a list so we need to remove this element as for it to not appear pending
                                    QPend.erase(i++);
                                    // If this was 'special_index's last appearance, mark it available for integration :-)
                                    // this is only valid if this so-called 'special_index' is different from the current one
                                    if (special_index != ix) {
                                        was_last_appearance = true;
                                        for (auto j = QPend.begin(); j != QPend.end(); ++j) {
                                            if (*j == special_index) {
                                                was_last_appearance = false;
                                            }
                                        }
                                        if (was_last_appearance) {
                                            PRINTF_DBG("We have effectively recieved the answer to all the sent msgs\n");
                                            std::cout << std::flush;
#pragma critical
{
                                            REF.p_READY_FOR_INTEGRATION->push(ix);
}
                                            ++total_processed;
                                        }
                                    } else {
                                        ++sent_locally;
                                    }
                                    std::cout << std::flush; // DEBUG :O
                                } else {
                                    ++i;
                                }
                            }
                            if (waiting){
                                answer_messages<0, 1, BATCH>(REF, O.MY_THREAD_n); // only one lap and no delay, i.e. TIMETOL = 1, DT = 0
                            }
                        }
                        PRINTF_DBG("DISPATCHED CORRECTLY\n");
                        std::cout << std::flush;
                    }
                }
            }
//            std::cout << "e" << std::endl;
//            std::cout << std::flush; // DEBUGGING
            // If we got answers to all sent messages then flag it globally as available for integration
            if (sent_locally == tot_locals){
                PRINTF_DBG("We have effectively recieved the answer to all the sent msgs\n");
                std::cout << std::flush;
#pragma critical
{
                REF.p_READY_FOR_INTEGRATION->push(ix);
}
                ++total_processed;
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

//            printf("Reporting so far having requested %d nodes... there are in total %d nodes\n",
//                   total_processed, atomic_helper);
//            std::cout << std::flush;

            ix_update = false;
        } else {
            tot_locals  = 0;
            sent_locally = 0;
            total_processed = 0;
            // this is what happens if there was no index update haha ;-)
            answer_messages<DT, TIMETOL, BATCH>(REF, O.MY_THREAD_n);
        }
//        std::cout << "g" << std::endl;
//        std::cout << std::flush; // DEBUGGING
        // (Re-)Attempt to grab a new index :K
#pragma omp critical
{
            if (!REF.p_CHECKED->empty()) {
                ix = REF.p_CHECKED->front();
                REF.p_CHECKED->pop();
                ix_update = true;
                PRINTF_DBG("Obtained a new ix\n");
                std::cout << std::flush;
            }
}
//        std::cout << "h" << std::endl;
//        std::cout << std::flush; // DEBUGGING
        // Add to the number of totals processed the ones from previous lap :-)
        if (total_processed != 0) {
#pragma omp atomic update
            *(REF.p_TOT) += total_processed;
            PRINTF_DBG("Increased the TOT from global pool\n");
            std::cout << std::flush;
        }
//        std::cout << "j" << std::endl;
//        std::cout << std::flush; // DEBUGGING

        // if there was no available index, check if it is still worth looping.
        if (!ix_update) {
#pragma omp atomic read
            atomic_helper = *(REF.p_TOT);
            globalstatus = (atomic_helper < NNodes);
        }
//        std::cout << "k" << std::endl;
//        std::cout << std::flush; // DEBUGGING
    }
    std::cout << "l" << std::endl;
    std::cout << std::flush; // DEBUGGING
    // One way or another being here means being done with the previous section :-)
    // If there are still pending messages then we should wait for them :-)
    total_processed = 0;
    waiting = (QPend.size()>0);
    while (waiting) {
        PRINTF_DBG("We are waiting before exiting, i.e. CLEANING! :-(\n");
        std::cout << std::flush;
        //mssleep(100);
        auto i = QPend.begin();
        while (i != QPend.end()) {
            MPI_Request_get_status(requests_send[*i], &sstatus, MPI_STATUS_IGNORE);
            MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
            if ((sstatus==1) && (rstatus==1)) {
                PRINTF_DBG("Recieved a response to our request: val %f\n", vval[*i]);
                std::cout << std::flush;
                std::get<0>(results[*i]) = vval[*i];
                special_index = ixlist[*i];
                (*REF.p_IntHelper)[special_index].ResultsPendProcess.push_back(results[*i]);
                QAvailable.push(*i);
                QPend.erase(i++);
                waiting = false;
                // Now reinitializing the requests won't be necessary :-)
//                requests_send[*i] = MPI_Request();
//                requests_recv[*i] = MPI_Request();
                was_last_appearance = true;
                for (auto j = QPend.begin(); j != QPend.end(); ++j) {
                    if (*j == special_index) {
                        was_last_appearance = false;
                    }
                }
                if (was_last_appearance) {
                    PRINTF_DBG("[LAST REMANENT] We have effectively recieved the answer to all the sent msgs\n");
                    std::cout << std::flush;
#pragma critical
                    {
                        REF.p_READY_FOR_INTEGRATION->push(ix);
                    }
                    ++total_processed;
                }
            } else {
                ++i;
            }
        }
    }
    // Update all the totals with the newfound information
#pragma omp atomic update
    *(REF.p_TOT) += total_processed;
} // end of function
















void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
