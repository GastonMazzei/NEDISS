//
// Created by m4zz31 on 5/11/21.
//

#ifndef CPPPROJCT_COMMUNICATIONFUNCTIONS_H
#define CPPPROJCT_COMMUNICATIONFUNCTIONS_H

#include <boost/serialization/string.hpp>
#include "../GraphClasses/GeneralGraph.h"
#include "../GraphClasses/GraphFunctions.h"
#include "../Utils/HelperClasses.h"



void GetAllMsgs(int NNodes, CommunicationHelper &H, Graph &g, ParallelHelper &P, IntegrationHelper &I, std::queue<long> &C);



void GetOneMsg(int ix,
               CommunicationHelper &H,
               Graph &g,
               ParallelHelper &P,
               IntegrationHelper &I,
               std::queue<long> &C);

void ask_for_node(int owner, VD &v, CommunicationHelper &H, int ix, Graph &g);

void ask_for_node_and_vertex(int owner, VD &v, ED &e, CommunicationHelper &H, int ix, Graph &g);

#endif //CPPPROJCT_COMMUNICATIONFUNCTIONS_H
