//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_GRAPHFUNCTIONS_H
#define CPPPROJCT_GRAPHFUNCTIONS_H
#include "../macros/macros.h"
#include "../Communication/CommunicationFunctions.h"
#include "../Solvers/GeneralSolver.h"
#include "../Utils/adequate_synchronization.h"
#include "../Utils/memory_management.h"
#include "../Utils/msleep.h"
#include "../Utils/error.h"
#include "GeneralGraph.h"
#include "../Utils/HelperClasses.h"
#include <mpi.h>
#include <set>
#include <boost/mpi/environment.hpp>
#include "../Communication/CommunicationFunctions.h"


template<typename DIFFEQ, typename SOLVER>
void contribute_to_integration(ReferenceContainer &REF, GeneralSolver<DIFFEQ,SOLVER> &solver, int MYPROC){

    PRINTF_DBG("Starting the contribution to int\n");std::cout<< std::flush;
    // Define timeout utilities
    const int DT= 5;
    int totlaps=0;
    // Build a graph vertex iterator
    auto vs = vertices(*REF.p_g);
    auto start = vs.first;
    auto end = vs.second;
    auto v = start;

    // Define relevant variables
    bool keepGoing = true;
    bool is_in = false;
    long * i;
    int remaining;
    long ix;
    unsigned long uix, currentuix;
    bool wait = false;
    std::vector<InfoVecElem> temp;

#pragma omp critical
    {
        remaining = *(REF.p_PENDING_INT);
    }


    // Main loop
    while (keepGoing) {
        // Fetch the data
#pragma omp critical
        {
            if (!REF.p_READY_FOR_INTEGRATION->second.empty()){
                ix = REF.p_READY_FOR_INTEGRATION->first.front();
                uix = REF.p_READY_FOR_INTEGRATION->second.front();
                REF.p_READY_FOR_INTEGRATION->first.pop();
                REF.p_READY_FOR_INTEGRATION->second.pop();
            } else {
                wait = true;
            }
        }

        if (wait) {
            mssleep(DT);
            PRINTF_DBG("(NOTHING TO INT)\n");std::cout<< std::flush;
        } else {
            ++totlaps;
            PRINTF_DBG("starting to locate this uix: %lu\n",uix);
            v = vertices(*REF.p_g).first;
            is_in =  false;
            while (!is_in) {
                if (get(get(boost::vertex_owner, *(REF.p_g)), *v) == MYPROC) {
                    currentuix = (unsigned long) get(get(boost::vertex_index, *(REF.p_g)), *v);
                    if (currentuix == uix) {
                        is_in = true;
                    } else {
                        v++;
                    }
                } else ++v;
                if (v == end) {
                    std::cout << "[CRITICAL] Failed to found ix!!!" << std::endl;
                };
            }
            // Join the data as required by the integrator
            temp = std::vector<InfoVecElem>();
            temp.resize((*REF.p_IntHelper)[ix].ResultsPendProcess.size());
            int _it = 0;
            auto it = (*REF.p_IntHelper)[ix].ResultsPendProcess.begin();
            auto end = (*REF.p_IntHelper)[ix].ResultsPendProcess.end();
            int oldsize = (*REF.p_IntHelper)[ix].ResultsPendProcess.size();
            std::set<unsigned long> seen;
            while (it != end){
                unsigned long elemix = std::get<2>(*it);
                const bool is_in = seen.find(elemix) != seen.end(); //CHECKPOINT
                PRINTF_DBG("Std output uix: %lu ix: %ld proc: %d. vals are %f , %f , %lu\n",
                       uix, ix, MYPROC,
                       std::get<0>(*it), std::get<1>(*it), std::get<2>(*it));
                std::cout<<std::flush;
                if (!is_in) {
                    temp[_it] = *it;
                    _it++;
                    seen.insert(elemix);
                } else {
                }
                (*REF.p_IntHelper)[ix].ResultsPendProcess.erase(it++);
            }
            temp.resize(_it);
            std::sort(temp.begin(),
                      temp.end(),
                      [](InfoVecElem &t1, InfoVecElem &t2) {
                          return std::get<2>(t1) > std::get<2>(t2);
                      }
            );
            PRINTF_DBG("Finished sorting! :-)\n");
            std::vector<double> neighborValues;
            std::vector<double> edgeValues;
            neighborValues.resize(temp.size());
            edgeValues.resize(temp.size());
            for (int j=0; j<temp.size(); ++j){
                neighborValues[j] = std::get<0>(temp[j]);
                edgeValues[j] = std::get<1>(temp[j]);
            }
            double result = -420.;
            solver.evolve(
                            (*REF.p_g)[*v].value,
                            (*REF.p_g)[*v].params,
                            neighborValues,
                            edgeValues,
                            result
                          );
            PRINTF_DBG("Integrating: temporal result was: %f", result);
            (*REF.p_g)[*v].temporal_register = result;
            PRINTF_DBG("Integrating: now it is: %f", (*REF.p_g)[*v].temporal_register);
        }
        wait = false;
#pragma omp critical
{
        remaining = REF.p_READY_FOR_INTEGRATION->first.size();
}

        PRINTF_DBG("There are %d remaining vals to integrate...\n", remaining);std::cout<<std::flush;

        // Recompute the status of 'KeepGoing'
        keepGoing = (remaining != 0);
    }
    PRINTF_DBG("In total this thread of PROC %d performed %d laps!\n",MYPROC, totlaps);std::cout<<std::flush;
}


template<int DT, int TIMETOL, int BATCH>
void answer_messages(ReferenceContainer &REF, int MYTHR);


template<int DT, int TIMETOL, int BATCH>
void answer_messages_edges(ReferenceContainer &REF,int MYTHR);

template<int DT, int TIMETOL, int BATCH>
void answer_field_requests(ReferenceContainer &REF,int MYTHR, int Nfield);

// Three new ones pend implementing :-)
template<int DT, int TIMETOL, int BATCH>
void perform_field_requests(ReferenceContainer &REF,int MYTHR, int fieldOrder);
template<typename DIFFEQ, typename SOLVER, int BATCH>
void contribute_to_higher_integration_without_communication(ReferenceContainer &REF, GeneralSolver<DIFFEQ,SOLVER> &solver, int fieldNum){};
template<typename DIFFEQ, typename SOLVER, int BATCH>
void finalize_integration(ReferenceContainer &REF, GeneralSolver<DIFFEQ,SOLVER> &solver){};

void sendReqForTest(int MYPROC, int i);

template <int DT, int TIMETOL, int BATCH>
void perform_requests(int NNodes,
                      ReferenceContainer REF,
                      unsigned long N,
                      OpenMPHelper &O);

// processors will iterate their "N=num_vertex(g)" nodes,
// and in each iteration they will either
//          (A)  (volume term)
//               iterate through their "NLocal=rank" locally-available neighbors
//               uploading the NLocal results into a container that can be used
//               to compute the central node's value using the solver evolver :-)
//          (B)  (surface term)
//               iterate through their "Nlocal<rank" locally-available neighbors,
//               while collecting "Nnonlocal" (what boost::graph::distributed call)
//               'vertex descriptors', 'edge descriptors', and the vertex's owner
//               in order to get involved in the following step which is requiring
//               a copy of this DynamicNode from the node's owner. Also, a new
//              iteration is performed over inward edges in order to collect the same
//              from them.
//
// It is in this context that we require keeping track of the surface (B) and volume (A) term's
// information using "N" (nodes in the process) containers that store a parallel cell object.

void register_to_value(Graph &g);

void destroyRequestWithoutCounter(MPI_Request &R);
void freeRequestWithoutCounter(MPI_Request &R);

// Fully declared template
template<typename DIFFEQ, typename SOLVER, int BATCH> // e.g. DIFFEQ = NoiselessKuramoto, SOLVER = EulerSolver<NoiselessKuramoto>
void single_evolution(Graph &g,
                      GeneralSolver<DIFFEQ,SOLVER> &solver,
                      CommunicationHelper &ComHelper,
                      ParallelHelper &ParHelper,
                      IntegrationHelper &IntHelper,
                      MappingHelper &MapHelper,
                      LayeredSolverHelper &LayHelper,
                      unsigned long N_total_nodes){

    // Roughly we are:
    // (1) gathering info from neighbors in parallel using MPI calls to other procs.
    // (2) solving the respective differential equation
    // (3) calling the proc synchronization

    unsigned long NVtot = boost::num_vertices(g), NT;
    int PENDING_INT = NVtot;


    std::pair<std::queue<long>, std::queue<unsigned long>> CHECKED;
    std::pair<std::queue<long>, std::queue<unsigned long>> READY_FOR_INTEGRATION;

    NT = ComHelper.NUM_THREADS;
    auto vs = vertices(g);
    double temporalResult;
    bool is_unclaimed = true;
    int TOT = 1;
    int active_responders=0, finalized_responders=0;
    bool keep_responding = true;
    ReferenceContainer REF(ParHelper,
                           ComHelper,
                           g,
                           CHECKED,
                           READY_FOR_INTEGRATION,
                           IntHelper,
                           TOT,
                           PENDING_INT,
                           MapHelper,
                           LayHelper,
                           keep_responding);

    const int MAX_SUBTHR = 1;
    const int TIMETOL = 20;
    const int DT = 1;

    bool isLayerBuilt = LayHelper.built;
    int request_performers=0;
    int request_performers_ended=0;


// The solver must be copied by each thread, as the 'FlowSpecs' attribute must act non-atomically.
#pragma omp parallel firstprivate(NVtot, vs, NT, N_total_nodes, REF, MAX_SUBTHR, TIMETOL, DT, solver, isLayerBuilt) // Node iterators have random access ;-)
    {
        bool am_i_first;
        int SplitCoef = omp_get_num_threads()/2; // 3/2
        if (SplitCoef < 2){SplitCoef = 2;}
        OpenMPHelper OmpHelper(NVtot, SplitCoef);
        if (OmpHelper.MY_THREAD_n < SplitCoef){
            if (OmpHelper.MY_THREAD_n % (MAX_SUBTHR + 1) == 0) {
                bool atomic_bool;
#pragma omp atomic read
                atomic_bool = keep_responding;
#pragma omp atomic update
                ++active_responders;
                    while (atomic_bool) { // as long as we keep processing our own,
                    //                    we mantain at least one dispatcher alive :-)
                        answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                        answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                        if (solver.requires_communication){
                            answer_field_requests<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n, 1);
                        }
#pragma omp atomic read
                        atomic_bool = keep_responding;
                        mssleep(DT);
                    }
#pragma omp atomic update
                ++finalized_responders;

                if (!atomic_bool) PRINTF_DBG("OVER!  :O\n");
                PRINTF_DBG(" I am thread %d (MESSAGE ANSWERER) and I have finished doing my job ;-)\n",omp_get_thread_num());

            } else {
#pragma omp atomic update
                request_performers++;
                perform_requests<DT, TIMETOL, BATCH>(NVtot, REF, N_total_nodes, OmpHelper);
#pragma omp atomic update
                request_performers_ended++;
            }
        } else {
            unsigned long NLocals, NInedges, M, rank, NOwned;
            long i=-1;
            unsigned long ui =0;

            bool ready4int = false;
            i += OmpHelper.MY_OFFSET_n;
            for (auto v = vs.first + OmpHelper.MY_OFFSET_n;
                 v != vs.first + OmpHelper.MY_OFFSET_n + OmpHelper.MY_LENGTH_n; ++v) {
                ++i;

                // Capture the central node's values
                PRINTF_DBG("Accesing values 1\n");std::cout<<std::flush;
                IntHelper[i].centralValue = g[*v].value;
                PRINTF_DBG("Accesing values 2\n");std::cout<<std::flush;
                IntHelper[i].centralParams = g[*v].params;
                IntHelper[i].build(g, *v, MapHelper, NOwned, rank, NLocals, M);
                if (NOwned == rank){
                    ready4int = true;
                } else {
                    ready4int = false;
                };
                if (get(MapHelper.NodeOwner,*v) != REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n]){
                    // we are treating a node which we dont own, please ignore it
                    PRINTF_DBG("\n\n\nIGNORING NODE THAT IS NOT OURS!\n\n\n\n");std::cout<<std::flush;
                    //++v;
                    exit(1);
                } else {
                    ui = ((unsigned long) get(get(boost::vertex_index, g), *v));
                }
                if (isLayerBuilt) LayHelper.buildForRank((long) ui, (long) rank);
                auto neighbors = boost::adjacent_vertices(*v, g);
                auto in_edges = boost::in_edges(*v, g);
                unsigned long MYPROCN = REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n];
#pragma omp parallel firstprivate(i, v, M, neighbors, NLocals, NInedges, NOwned, N_total_nodes, MYPROCN)
                {
                    OpenMPHelper OmpHelperN(NLocals, 0);
                    for (auto n = neighbors.first + OmpHelperN.MY_OFFSET_n;
                         n != neighbors.first + OmpHelperN.MY_OFFSET_n + OmpHelperN.MY_LENGTH_n; ++n) {
                        auto local_e = edge(*v, *n, g).first;
                        auto local_v = *n;
                //if (get(MapHelper.Local, *n) == 1) { // Case (1): locally available elements ;-)
                if (REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n] == get(get(boost::vertex_owner, g), local_v)) {
                        if (edge(*v, *n, g).second == 1) {
                                PRINTF_DBG("Accesing values 3: are we the owner? Nowner: %d, we: %d, Vowner %d, Eowner: %d (<-temp placeholder) MapHelper.Local applied to neighbor yields: %d\n",
						get(get(boost::vertex_owner, g), local_v),
						process_id(g.process_group()),
						get(get(boost::vertex_owner, g),*v),
						0,
						get(MapHelper.Local, *n));std::cout<<std::flush;

#pragma omp critical         // Just store the results directly in the results.
                                // we can probably afford being critical as there is
                                // MPI communication overhead in other threads.
                                {
                                    IntHelper[i].ResultsPendProcess.emplace_back(g[*n].value,
                                                                                 g[edge(*v, *n,
                                                                                        g).first].value, // index to UID :-)
                                    (unsigned long) ((unsigned long) get(get(boost::vertex_index, g), *n) + N_total_nodes * ((unsigned long) MYPROCN) ) );
                                }
                            } else { error_report("Push back mechanism for local nodes has failed"); };
                        } else { // Case (2.A): We 'see' this neighbor because it is connected
                            // via an edge that we own. Get the edge and record who is the other node's owner
                            PRINTF_DBG("Accesing values 4\n");std::cout<<std::flush;

#pragma omp critical
                            {
                                ParHelper.data[i].MissingA[OmpHelperN.MY_THREAD_n].emplace_back(g[local_e].value,
                                                                                                (int) get(MapHelper.NodeOwner,
                                                                                                    *n),
                                                                                                (int) get(get(boost::vertex_index,
                                                                                                        g), *n));
                            }
                        }
                    }
                    // Finally please contemplate the 3rd case:
                    // we have in-edges that are not available
                    // locally so we have neighbors that cant
                    // be indexed by 'neighbors'
                    OpenMPHelper OmpHelperE(M, 0, OmpHelperN.N_THREADS_n, OmpHelperN.MY_THREAD_n);
                    int j = 0;
                    long central_ix;
#pragma omp critical
                    {
                        central_ix = get(get(boost::vertex_index,g), *v);
                    }
                    for (auto e = in_edges.first; e != in_edges.second; ++e) {
                        if ((j >= OmpHelperE.MY_OFFSET_n) && (j < OmpHelperE.MY_OFFSET_n + OmpHelperE.MY_LENGTH_n)) {
                            auto local_e = *e;
                            auto local_v = boost::source(*e, g);
#pragma omp critical
                            {
                                ParHelper.data[i].MissingB[OmpHelperE.MY_THREAD_n].emplace_back(
                                        (double) central_ix,
                                        (int) get(MapHelper.NodeOwner,
                                            local_v),
                                        (int) get(get(boost::vertex_index,
                                                g), local_v));
                            }
                        }
                        ++j;
                    }
                } // end of the nested parallelism! :-)


                if (ready4int){
#pragma omp critical
                    {
                        READY_FOR_INTEGRATION.first.push(i);
                        READY_FOR_INTEGRATION.second.push(ui);
                    }
#pragma omp atomic update
                    ++TOT;
                    PRINTF_DBG("Increased by one the number of total vertex\n");std::cout<<std::flush;
                } else {
#pragma omp critical
                    {
                        CHECKED.first.push(i); // Adding the index to the list of checked indexes ;-)
                        CHECKED.second.push(ui); // Adding the index to the list of checked indexes ;-)
                    }
                    PRINTF_DBG("Added ix %d to CHECKED\n", i);std::cout<<std::flush;
                }
            }
            PRINTF_DBG(" I am thread %d (FOR WORKER) and I have finished doing my job ;-)\n",omp_get_thread_num());
            PRINTF_DBG("This thread is a for worker that has ended :-)\n");
        }

        PRINTF_DBG("about to enter 'omp critical unclaimed'\n"); std::cout << std::flush;
#pragma omp critical
        {
            am_i_first = is_unclaimed;
            if (is_unclaimed){
                is_unclaimed = false;
                PRINTF_DBG("Claimed by thread %ld of processor %d\n", OmpHelper.MY_THREAD_n, ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n]);
            }
        }

        if (am_i_first) {
            // Prepare vars
            int atomical_int;
            bool are_we_over = false;

            // While our Process is not over, I am the official responder :-)
#pragma omp atomic read
            atomical_int = TOT;
            PRINTF_DBG("Already checking if we are over, proc %d\n", ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n]);
            are_we_over = (atomical_int >= NVtot); // TODO: bug. it should be exactly equal (looks at perform_requests)
            int notreadyyet=0;

            // FIRST are we over: checking that all nodes have been processed
            while (!are_we_over){
                // Re-check if we are over
#pragma omp atomic read
                atomical_int = TOT;
                are_we_over = (atomical_int >= NVtot); // TODO: same here ;-).
                ++notreadyyet;
                mssleep(DT);
                PRINTF_DBG("not ready yet!\n");
            }

            // SECOND are we over: checking that all perform requesters have exited
            int aux1,aux2;
#pragma omp atomic read
            aux1 = request_performers;
#pragma omp atomic read
            aux2 = request_performers_ended;
            are_we_over = (aux1 == aux2);
            while (!are_we_over){
#pragma omp atomic read
                aux2 = request_performers_ended;
                are_we_over = (aux1 == aux2);
                mssleep(DT);
            }

            PRINTF_DBG("\n\n\n\n\n\n\WE WERE OVERRR NVtot and TOT are: %d & %d\n\n\n\n", NVtot, TOT);
            PRINTF_DBG("Before being ready, we waited for %d laps!\n", notreadyyet);
            std::cout << std::flush;


            int worldsize = ComHelper.WORLD_SIZE[OmpHelper.MY_THREAD_n];
            int worldrank = ComHelper.WORLD_RANK[OmpHelper.MY_THREAD_n];

            PRINTF_DBG("\nabout to 1st barrier says P:%d\n", worldrank);std::cout<<std::flush;
            MPI_Barrier(MPI_COMM_WORLD);

            PRINTF_DBG("The (FIRST) cherry of the cake ;-)\n");std::cout<<std::flush;
#pragma omp atomic write
            keep_responding = false;

            bool readme;
#pragma omp atomic read
            readme = keep_responding;
            PRINTF_DBG("keep_responding was now set as %d\n",readme);std::cout<<std::flush;

            int auxy1, auxy2;
#pragma omp atomic read
            auxy1 = active_responders;
#pragma omp atomic read
            auxy2 = finalized_responders;
            int jh=0,jmax=1000;
            while (auxy1 > auxy2){
                mssleep(DT);
#pragma omp atomic read
                auxy2 = finalized_responders;
#pragma omp atomic read
                auxy1 = active_responders;
                PRINTF_DBG("active (total) are %d and finalized are %d\n", auxy1, auxy2);
                jh++;
                if (jh>jmax) PRINTF_DBG("active (total) are %d and finalized are %d\n", auxy1, auxy2);
            }
            PRINTF_DBG("\nabout to 2nd barrier says P:%d\n", worldrank);std::cout<<std::flush;
            MPI_Barrier(MPI_COMM_WORLD);
            PRINTF_DBG("\nEnded barrier says P:%d\n", worldrank);std::cout<<std::flush;

            PRINTF_DBG("The (SECOND) cherry of the cake (-;\n");std::cout<<std::flush;

        } else if (OmpHelper.MY_THREAD_n % 2 == 0) {
            bool atomic_bool;
#pragma omp atomic read
            atomic_bool = keep_responding;
            bool later_mark_finalized = false;
            if (atomic_bool){
#pragma omp atomic update
                active_responders ++;
                later_mark_finalized = true;
                }
            while (atomic_bool){
                answer_messages<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                answer_messages_edges<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n);
                if (solver.requires_communication) {
                    answer_field_requests<DT, TIMETOL, BATCH>(REF, OmpHelper.MY_THREAD_n, 1);
                }
                // 2) Re-check if we are over
#pragma omp atomic read
                atomic_bool = keep_responding;
                PRINTF_DBG("so far I keep responding U.u cuz keep responding was %d\n", atomic_bool);
            }
            if (later_mark_finalized) {
#pragma omp atomic update
                finalized_responders++;
            }
        }

/// COULD THIS BE INSERTED INSIDE THE PREVIOUS PARALLEL REGION? todo: optimization upgrade.

        // Integration section: this is only exited when there are no more pending integration remaining :-)
        PRINTF_DBG("starting to contribute\n");std::cout<<std::flush;
        int auxv1,auxv2;
#pragma omp atomic read
        auxv1 = request_performers;
#pragma omp atomic read
        auxv2 = request_performers_ended;
        bool keep_integrating = (auxv1 > auxv2);
        while (keep_integrating) {
            contribute_to_integration<DIFFEQ, SOLVER>(REF, solver,
                                      REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n]);
#pragma omp atomic read
            auxv1 = request_performers;
#pragma omp atomic read
            auxv2 = request_performers_ended;
            keep_integrating = (auxv1 > auxv2);
            mssleep(DT);
        }
        contribute_to_integration<DIFFEQ, SOLVER>(REF, solver,
                                      REF.p_ComHelper->WORLD_RANK[OmpHelper.MY_THREAD_n]);
        PRINTF_DBG("ending to contribute\n");std::cout<<std::flush;


} // end of the parallel construct

printf("Exited the parallel construct! yay! solver deg is %d\n", solver.deg);

    for (int i=2; i < solver.deg ; ++i){
        printf("entering the for, i is %d\n",i);
        if (solver.requires_communication){
	        bool keep_responding = true;
            int Ncapturers = 0;
            int NFinalizedCapturers = 0;
            int Nresponders = 0;
            int NFinalizedResponders = 0;
            long pending = (long) NVtot;
            std::queue<long> CAPTURED; 
#pragma omp parallel
{
            if ((omp_get_thread_num() % 2 == 1) && (omp_get_thread_num != 0)){
#pragma omp atomic update
                Ncapturers++;
                long locallyPending;
#pragma omp atomic read
                locallyPending = pending;
                while (locallyPending > 0){
                    perform_field_requests<DT, TIMETOL, BATCH>(REF,
                     REF.p_ComHelper->WORLD_RANK[omp_get_thread_num()],
                     i);
#pragma omp atomic read
                    locallyPending = pending;
                    mssleep(DT);
                    // -- only -- for -- debugging -- haha
#pragma omp atomic update
                    --pending;
                    printf("locallyPending' is now: %l\n",locallyPending);
                    //----end--of--debugging--section--hahaa
		         }
#pragma omp atomic update
		        NFinalizedCapturers++;
	         } else if ((omp_get_thread_num() % 2 == 0) && (omp_get_thread_num != 0)){
#pragma omp atomic update
		        Nresponders++;
                bool l_keep_responding = true;
		        while (l_keep_responding){
                    answer_field_requests<DT, TIMETOL, BATCH>(REF, (int) omp_get_thread_num(), i);
#pragma omp atomic read
                    l_keep_responding = keep_responding;
                    mssleep(DT);
		        }
#pragma omp atomic update
		        NFinalizedResponders++;
            } else {
                int v_Ncapturers = 0;
                int v_NFinalizedCapturers = 0;
                int v_Nresponders = 0;
                int v_NFinalizedResponders = 0;
                long locallyPending;
#pragma omp atomic read
                v_Ncapturers = Ncapturers;
                while (v_Ncapturers == 0){
#pragma omp atomic read
                    v_Ncapturers = Ncapturers;
                }
#pragma omp atomic read
                v_NFinalizedCapturers = NFinalizedCapturers;
                while (v_Ncapturers > v_NFinalizedCapturers){
#pragma omp atomic read
                    v_Ncapturers = Ncapturers;
#pragma omp atomic read
                    v_NFinalizedCapturers = NFinalizedCapturers;
                    mssleep(DT);
                    printf("They are not yet equal: %d vs %d\n", (int) v_Ncapturers, (int) v_NFinalizedCapturers);
                }
                printf("about to implement the first integration barrier\n");
                MPI_Barrier(MPI_COMM_WORLD);
#pragma omp atomic write
                keep_responding = false;
#pragma omp atomic read
                v_Nresponders = Nresponders;
                while (v_Nresponders == 0){
#pragma omp atomic read
                    v_Nresponders = Nresponders;
                }
#pragma omp atomic read
                v_NFinalizedResponders = NFinalizedResponders;
                while (v_Nresponders > v_NFinalizedResponders){
#pragma omp atomic read
                    v_Nresponders = Nresponders;
#pragma omp atomic read
                    v_NFinalizedResponders = NFinalizedResponders;
                    mssleep(DT);
                }
                MPI_Barrier(MPI_COMM_WORLD);
	        }
}
	// compute the next required step
        } else {
		    contribute_to_higher_integration_without_communication<DIFFEQ, SOLVER, BATCH>(REF,
                                                                    solver,
				                                                    i);
	    }
    }
    // Join all the integration terms
    printf("Reached the final integration :-)\n");
    finalize_integration<DIFFEQ, SOLVER, BATCH>(REF, solver);
    


    PRINTF_DBG("starting to swap register\n");std::cout<<std::flush;
    // Swap temporal and main registers
    register_to_value(g);

    PRINTF_DBG("exited");
    PRINTF_DBG("About to synchronize");std::cout<<std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    PRINTF_DBG("About to increase time by h\n");
    solver.EvolveTime();
    PRINTF_DBG("Done");
    printf("-^-^-H-E-A-R-T-B-E-A-T-^-^-");
    PRINTF_DBG("\n\n\n\\n\n\n\n\n\n");std::cout<<std::flush;
}





#endif //CPPPROJCT_GRAPHFUNCTIONS_H
