//
// Created by m4zz31 on 5/11/21.
//


#include "CommunicationFunctions.h"


void build_answer_edges(double * answer, ReferenceContainer &REF, double * ix, int owner, int &MyNProc){
    bool found = false;
    auto vs = vertices(*REF.p_g);
    auto v = vs.first;
    while ((v != vs.second) && (!found)){
        if (get(REF.p_MapHelper->NodeOwner,*v)==MyNProc) {

            // If I am the owner, test if this node has the index I'm looking for
            if (get(get(boost::vertex_index, *(REF.p_g)), *v) == (int) *(ix+1)) {

                // If this is the node this code was aimed to, iterate through neighbrors
                auto neighbors = boost::adjacent_vertices(*v, *REF.p_g);
                auto n = neighbors.first;
                while ((n != neighbors.second) && (!found)) {

                    //Iterate through neighbors looking for the correct edge ;-)
                    if (get(REF.p_MapHelper->NodeOwner, *n) == owner) {

                        // If this neighbor is owned by the requester
                        if (get(get(boost::vertex_index, *(REF.p_g)), *n) == (int) *ix) {

                            // If its index is the one the requester refered to
                            found = true;
                            PRINTF_DBG("Found the requested node attached to the requested edge!\n");
                            auto e = edge(*v, *n, *(REF.p_g));
                            if (e.second != 1) {
                                printf("[FATAL] Requested edge did not exist!!!");
                                exit(1);
                            }
                            int TAG = OFFSET + (int) *ix;
                            *answer = (*(REF.p_g))[*v].value;
                            *(answer + 1) = (*(REF.p_g))[e.first].value;
                        }
                    }
                    ++n;
                }
            }
        }
        ++v;
    }
    if (!found) {
        printf("CRITICAL: node+edge not found!\n");
        std::cout << "We were asked for " <<  *ix << " and " << *(ix+1) << std::flush;
        std::cout << "owner and mynproc are " <<  owner << " and " << MyNProc << std::flush;
        PRINTF_DBG("[CRITICAL] THE NODE WAS NOT FOUND\n");
        std::cout << std::flush;
        exit(1);
    }
}


void build_answer_nodes(double &answer, ReferenceContainer &REF, double ix, int owner, int &MyNProc){
    auto vs = vertices(*REF.p_g);
    for (auto v = vs.first; v != vs.second; ++v){
        if (get(REF.p_MapHelper->NodeOwner,*v)==MyNProc) {
            if (get(get(boost::vertex_index, *(REF.p_g)), *v) == (int) ix){
                answer = (*REF.p_g)[*v].value;
                return;
            }
        }
    }
    PRINTF_DBG("[CRITICAL] THE NODE WAS NOT FOUND\n");std::cout<<std::flush;
    exit(1);
}

















