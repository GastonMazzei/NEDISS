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
            if (flagprobes[i] == 1) { // a message appears to be available! recieve it!

                // Print the first three elements in the asynchronous probe :-)
                PRINTF_DBG("|%d %d %d|\n",flagprobes[0] , flagprobes[1] , flagprobes[2]);std::cout<<std::flush;

                // Try to recieve it
                R_tot.push_back(MPI_Request());
                recv_nonblocking(i, R_tot[R_tot.size() - 1], ix, VERTEXVAL_REQUEST_FLAG);
                statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
                while (statusreq_status != 0) {
                    statusreq_status = MPI_Test(&R_tot[R_tot.size() - 1], &status_localreq, MPI_STATUS_IGNORE);
                }


                if (status_localreq == 1) { // If we were first to capture the message, proceed.
                    R_send.push_back(MPI_Request());
                    PRINTF_DBG("We effectively captured a vertex info request :-)\n");
                    irespond_value(REF, ix, i, R_send[R_send.size() - 1], MYPROC);
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
            probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
            while (probe_status!=0){
                probe_status = MPI_Iprobe(i, VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD, &flagprobes[i], MPI_STATUS_IGNORE);
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
            PRINTF_DBG("ENDing answer_messages. Caught (R_tot)=%d requests, Answered (R_send)=%d.\n",
                       R_tot.size(),R_send.size());

            int status_of_getstatus = 1;
            for (int i = 0; i < R_send.size(); ++i) {
                status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
                while (status_of_getstatus != 0) {
                    printf("Failing to test a request, its returning %d\n", status_of_getstatus);
                    status_of_getstatus = MPI_Test(&R_send[i], &status_localreq, MPI_STATUS_IGNORE);
                }

                if (status_localreq == 0) {
                    PRINTF_DBG("We are currently waiting for a send to be completed\n");
                    MPI_Wait(&R_send[i], MPI_STATUS_IGNORE);
                    PRINTF_DBG("Successfully waited the request to be completed\n");
                }
            }
        }

        std::vector<MPI_Request> R_tot;
        std::vector<MPI_Request> R_send;
        mssleep(DT);
        ++ticks;
    }
    printf("---answer_request heartbeat---\n");std::cout<<std::flush;
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
    //MPI_Request requests_send[BATCH]; THIS HAS BEEN DISABLED AFTER THE LAST MIGRATION
    MPI_Request requests_recv[BATCH];
    InfoVecElem results[BATCH];
    std::queue<int> QAvailable;
    std::list<int> QPend;
    for (int i = 0; i < BATCH; ++i){
        QAvailable.push(i);
    }
    double vval[BATCH], their_vix[BATCH];
    int ixlist[BATCH]; // this is our new helper in order to process various indexes asynchronously (=^-|)-|-/
    int owner[BATCH];
    int retStatus=1;
    int status_rstatus=1;
    int counter=0,MAX_TRIES=20; // #HYPERPARAMS #HYPERPARAMETERS


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

            // race-condition unfriendly? Only extensive testing will help us decide :^|
            for (auto _it= itBeg; _it != itEnd; ++_it) {

                PRINTF_DBG("About to ask for the process of local elem ix: %d\n",ix);

                auto &thread = *_it;
                tot_locals += itEnd - itBeg;


                for (auto it =  thread.begin();it != thread.end(); it++){

                    // retrieve data
                    owner[QAvailable.front()] = std::get<1>(*it);
                    their_vix[QAvailable.front()] = std::get<2>(*it);

                    // Debug Station
                    PRINTF_DBG("Sending and recieving (nonblocking): asking %d for vertex w index: %f \n",
                           owner[QAvailable.front()], their_vix[QAvailable.front()]);
                    std::cout  << std::flush;



                    // ***NEW WAY*** A
                    // Send in a way that guarantees us that the message has been sent :-)
                    // This should slow the script down but will 'increase' the guarantees of messages arriving
                    while (retStatus != 0){
                        randi = gen(rng);
                        printf("About to perform MPI_Ssend w id %d\n", randi);
                        retStatus = MPI_Ssend(&their_vix[QAvailable.front()], 1, MPI_DOUBLE, owner[QAvailable.front()], VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD);
                        printf("Correctly performed MPI_Ssend w id %d\n", randi);
                        std::cout << std::flush;
                    }
                    retStatus = 1;
                    // ***OLD WAY*** A
//                    send_nonblocking(owner[QAvailable.front()],
//                                     requests_send[QAvailable.front()],
//                                     their_vix[QAvailable.front()], VERTEXVAL_REQUEST_FLAG);



                    // THIS IS CORRECT:

                    // Recieve the answer for the requested info, using tag = index
                    recv_nonblocking(owner[QAvailable.front()], // recieve the respone (tag=ix) for that  ix
                                     requests_recv[QAvailable.front()], // flag = (int) index :-)
                                     vval[QAvailable.front()], (int) their_vix[QAvailable.front()]);

                    // Store the result in our BATCH-sized temporal container:
                    // we are only lacking the exact vertex value which should be returned from the owner
                    results[QAvailable.front()] = std::make_tuple(0, // placeholder until we get the correct val
                                                                  std::get<0>(*it), // we are inaugurating this indexing model [:<)
                                                                  owner[QAvailable.front()] * N + their_vix[QAvailable.front()]);

                    QPend.push_back(QAvailable.front());
                    std::cout << std::flush;
                    ixlist[QAvailable.front()] = ix;
                    std::cout << std::flush;
                    QAvailable.pop();


                    // If QAvailable is empty, devote ourselves to asynchronously waiting for
                    // answers to arrive. This indirectly means that our QPend has reached length BATCH.
                    if (QAvailable.empty()) {
                        waiting = true;

                        // Debug Station
                        PRINTF_DBG("it IS FINALLY EMPTY\n");
                        std::cout << std::flush;


                        while (waiting) {

                            PRINTF_DBG("We are waiting to be responded :-( Qpends size is %d\n", QPend.size());
                            std::cout << std::flush;

                            // Iterate through the pending indexes to see if one has been answered
                            auto i = QPend.begin();
                            while (i != QPend.end()) {

                                rstatus = 0;
                                status_rstatus = 1;
                                counter=0;
                                bypass = false;
                                retStatus = 1;

                                // Waiting for the answer to arrive with a timeout.
                                while ((rstatus != 1) && (!bypass)){
                                    status_rstatus = MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                                    while (status_rstatus!=0){
                                        PRINTF_DBG("fetching the status_rstatus failed as it yielded %d:0\n", status_rstatus);std::cout<<std::flush;
                                        status_rstatus = MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                                    }
                                    ++counter;
                                    if ((counter>=MAX_TRIES) && (rstatus!=1)){ // we give one last chance to rstatus ;-D
                                        // the message probably didnt arrive, so we will resend it :-)
                                        bypass = true;
                                        while (retStatus != 0){
                                            retStatus = MPI_Ssend(&their_vix[*i], 1, MPI_DOUBLE, owner[*i], VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD);
                                        }
                                    }
                                }

                                if (rstatus == 1) {
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
                                    //MPI_Request_free(&requests_send[*i]);
                                    MPI_Request_free(&requests_recv[*i]);
                                    //requests_send[*i] = MPI_Request();
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
                                //answer_messages<0, 1, BATCH>(REF, O.MY_THREAD_n); // #HYPERPARAMS #HYPERPARAMETERS
                            }
                        }
                        PRINTF_DBG("DISPATCHED CORRECTLY\n");
                        std::cout << std::flush;
                    }
                }
            }

            if (sent_locally == tot_locals){
                PRINTF_DBG("We have effectively recieved the answer to all the sent msgs\n");
                std::cout << std::flush;
#pragma critical
{
                REF.p_READY_FOR_INTEGRATION->push(ix);
}
                ++total_processed;
            }

//          *******************THE SAME SHOULD BE DONE FOR EDGES PLEASE********************
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
//          *******************THE SAME SHOULD BE DONE FOR EDGES PLEASE********************



        } else {
            // *************THIS IS WHAT HAPPENS IF THERE WAS NO INDEX UPDATE****************
            tot_locals  = 0;
            sent_locally = 0;
            total_processed = 0;
            // this is what happens if there was no index update haha ;-)
            //answer_messages<DT, TIMETOL, BATCH>(REF, O.MY_THREAD_n); // #HYPERPARAMS #HYPERPARAMETERS
        }


//        ***SECOND SECTION: GRABBING NEW INDEXES AND CHECKING IF LOOPING STILL MAKES SENSE*******
        ix_update = false;
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


        // Add to the number of totals processed the ones from previous lap :-)
        if (total_processed != 0) {
#pragma omp atomic update
            *(REF.p_TOT) += total_processed;
            PRINTF_DBG("Increased the TOT from global pool\n");
            std::cout << std::flush;
        }


        // if there was no available index, check if it is still worth looping.
        if (!ix_update) {
#pragma omp atomic read
            atomic_helper = *(REF.p_TOT);
            globalstatus = (atomic_helper < NNodes);
        }

    }

    // ***************THIRD SECTION: END OF GLOBALSTATUS TRUE*****************************

    // Untidy Debug Station
//        std::cout << "k" << std::endl;
//        std::cout << std::flush; // DEBUGGING
    std::cout << "-3rd section.-\n" << std::endl;
    std::cout << std::flush; // DEBUGGING




    //          ********PROCESS REMAINING OPEN REQUESTS***************
    // One way or another being here means being done with the previous section :-)
    // If there are still pending messages then we should wait for them :-)
    //
    waiting = (QPend.size()>0);
    while (waiting) {

        // Iterate through the pending indexes to see if one has been answered
        auto i = QPend.begin();
        while (i != QPend.end()) {

            // Initialize vars :-)
            rstatus = 0;
            status_rstatus = 1;
            counter=0;
            bypass = false;
            retStatus = 1;

            // Waiting for the answer to arrive with a timeout.
            while ((rstatus != 1) && (!bypass)){
                status_rstatus = MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                while (status_rstatus!=0){
                    PRINTF_DBG("fetching the status_rstatus failed as it yielded %d:0\n", status_rstatus);std::cout<<std::flush;
                    status_rstatus = MPI_Request_get_status(requests_recv[*i], &rstatus, MPI_STATUS_IGNORE);
                }
                ++counter;
                if ((counter>=MAX_TRIES) && (rstatus!=1)){ // we give one last chance to rstatus ;-D
                    // the message probably didnt arrive, so we will resend it :-)
                    bypass = true;
                    while (retStatus != 0){
                        retStatus = MPI_Ssend(&their_vix[*i], 1, MPI_DOUBLE, owner[*i], VERTEXVAL_REQUEST_FLAG, MPI_COMM_WORLD);
                    }
                }
            }

            if (rstatus == 1) {
                // We keep track of the locally successfully sent and recieved
                // requests as to at last check if, for this index, we were able
                // to settle the issue and its ready for integration or not :-)
                std::get<0>(results[*i]) = vval[*i];
                special_index = ixlist[*i];
                (*REF.p_IntHelper)[special_index].ResultsPendProcess.push_back(results[*i]);

                // Mini debug station
                PRINTF_DBG("Recieved a response to our request: val %f for ix: %d\n", vval[*i], special_index);
                std::cout << std::flush;

                // Kill forever the MPI Request objects ;-)
                MPI_Request_free(&requests_recv[*i]);

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
                }
            } else {
                ++i;
            }
        }

        waiting = (QPend.size()>0);
    }
    PRINTF_DBG("DISPATCHED CORRECTLY THE REMAINING UNANSWERED REQUESTS :-)\n");
    std::cout << std::flush;

    // SUM TO THE GLOBAL TOTAL-OF-PROCESSED (i.e. "TOT")
#pragma omp atomic update
    *(REF.p_TOT) += total_processed;

    printf("Final termination of perform_requests :-)\n");std::cout<<std::flush;
} // end of function
















void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
