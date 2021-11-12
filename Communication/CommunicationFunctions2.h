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









    // Collect missing nodes:
    // Access to "MissingA" can be direct as there are not
    // race conditions :-)
    //
    // ITERATE OVER IX HERE
    //

// BEFORE EXITING PLEASE process the remaining cases :-) N :-)
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
    long current_ix;

    bool all_sent_locally = false;
    int tot_locals = 0;
    int sent_locally = 0;

    // insert here required initialization
    std::list<long> all_indexes_seen;
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
    std::pair<int, int> ixlist[BATCH]; // this is our new helper in order to process various indexes asynchronously (=^-|)-|-/ 
    int owner[BATCH];
    // end of the required initialization



#pragma omp atomic read
    atomic_helper = *(REF.p_TOT);
    globalstatus = (atomic_helper < NNodes);

    while (globalstatus) {
        // Perform the main task: ask for the required information!
            if (ix_update) {
                //mssleep(1000);
                //std::cout << "[DEBUG] about to call GetOneMsg" << std::endl;
                printf("About to ask for the process of local elem ix: %d\n",ix);

// BEGGINING OF THE MANUALLY INSERTED GetOneMsg core
// plz indent it :-) 
    auto itBeg = REF.p_ParHelper->data[ix].MissingA.begin();
    auto itEnd = REF.p_ParHelper->data[ix].MissingA.end();
    tot_locals  = 0;
    sent_locally = 0;
    for (auto _it= itBeg; _it != itEnd; ++_it) { // race-condition unfriendly?
        auto &thread = *_it;
	tot_locals += itEnd - itBeg;
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
	    // add to all_indexes[QAvailable.front()] the pair <QAvailable.front(), ix> 
	    // so that if we dont get the answer in this ix then we can do it in the future 
	    // setting 
	    // local_ix = all_indexes[*i].second ;-) 
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
			    ++sent_locally;
                            std::cout << std::flush;
                            std::get<0>(results[*i]) = vval[*i];
                            // No lock is needed: this index is entirely ours :-)i
			    // WARNING! this ix is not just 'ix' but the one from the std::pair<int,int>[BATCH]
			    //local_ix 
                            (*REF.p_IntHelper)[local_ix].ResultsPendProcess.push_back(results[*i]);
                            QAvailable.push(*i);
                            QPend.erase(i++);
                            waiting = false;
                            requests_send[*i] = MPI_Request();
                            requests_recv[*i] = MPI_Request();
			    // WARNING! check local_ix's status in all BATCH-1 elements in std::pair<int,int>[BATCH]
			    // if this index was the last one pending and its not the current, add to the 
			    // pend integration list and remove from the all_indexes std::list<long>
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

// END OF THE MANUALLY INSERTED GetOneMsg core
//
                PRINTF_DBG("Skipped sending messages :-)\n");
                ++total_processed; // this is O.K., I mean is kinda qualitative

                // Also communicate to other threads that there is one new item pending integration :o
		// but only do it if we have effectively recieved responses to all sent messages!
		// if (all_sent_locally) // and then set all_sent_locally to false ;-D
		if (sent_locally == tot_locals){
#pragma critical
                {
                    REF.p_READY_FOR_INTEGRATION->push(ix);
                }
		} else {
		// add to the list of indexes or potentially a set
		} 


                printf("Reporting so far having requested %d nodes... there are in total %d nodes\n",
                       total_processed, atomic_helper);
                ix_update = false;
            } else {
		    // this is what happens if there was no index update haha ;-)
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

// SOMEHOW NOW WE GOTTA REDEFINE GLOBAL STATUS SO THAT WE HAVE A CHECKER INDEPENDENT OF TOT, which we build only at the end.

#pragma omp atomic read // check if this loop still makes sense
            atomic_helper = *(REF.p_TOT);
            globalstatus = (atomic_helper < NNodes);
        }

// NOW ALLOW some sent messages to arrive
// NOTE: here we are using "ix", which at this point should not be defined... maybe not use 'ix' 
// but some other mechanism :-) 
//
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

    //add the list of indexes to the pend integration list :-).


// IMPORTANT!
// only after waiting for the entire arrival of all of our messages it is that we can do:
//  *(REF.p_TOT) += total_processed;

} // end of function















void ask_for_node(int owner, double &vvalue, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, double &vvalue, double &evalue, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
