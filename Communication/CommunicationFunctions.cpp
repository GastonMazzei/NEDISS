//
// Created by m4zz31 on 5/11/21.
//


#include "CommunicationFunctions.h"

#include "../Utils/global_standard_messages.h"

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
        std::cout << std::flush;
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





void FieldRequestObject::buildSendTag(int * data){
    sendTag = *data;
};

void FieldRequestObject::buildRecvTag(int * data){
    if (field == 0){
        recvTag = K1_REQUEST;
        //std::cout << "rk1 req tag is " << sendTag << std::endl;
        return;
    } else if (field == 1) {
        recvTag = K2_REQUEST;
        //std::cout << "rk2 req tag is " << sendTag << std::endl;
        return;
    } else if (field == 2) {
        recvTag = K3_REQUEST;
        //std::cout << "rk3 req tag is " << sendTag << std::endl;
        return;
    }
};


void FieldRequestObject::computeAnswer(ReferenceContainer &REF, int * buffer, double * answer){
    if (field == 0){
#pragma omp atomic read
        (*answer) = REF.p_LayHelper->data[*(buffer + 1)].RK1[0];
    } else if (field == 1) {
#pragma omp atomic read
        (*answer) = REF.p_LayHelper->data[*(buffer + 1)].RK2[0];
    } else if (field == 2) {
#pragma omp atomic read
        (*answer) = REF.p_LayHelper->data[*(buffer + 1)].RK3[0];
    }
};

void FieldRequestObject::computeReady(ReferenceContainer &REF, int * buffer, bool &isReady){
    if (field == 0){
#pragma omp atomic read
        isReady = REF.p_LayHelper->data[*(buffer + 1)].RK1_status;
        std::cout << "rk1 status returned " << isReady << std::endl;
        return;
    } else if (field == 1) {
#pragma omp atomic read
        isReady = REF.p_LayHelper->data[*(buffer + 1)].RK2_status;
        return;
    } else if (field == 2) {
#pragma omp atomic read
        isReady = REF.p_LayHelper->data[*(buffer + 1)].RK3_status;
        return;
    }
};

