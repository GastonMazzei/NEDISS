//
// Created by m4zz31 on 31/10/21.
//

#ifndef CPPPROJCT_ADEQUATE_SYNCHRONIZATION_H
#define CPPPROJCT_ADEQUATE_SYNCHRONIZATION_H

#include "../GraphClasses/GeneralGraph.h"

void adsync_message(std::string message, Graph &g);

void adsync_barrier();

void adsync_message_barrier(std::string detail, Graph &g);

void adsync_synchronization_barrier(std::string detail, Graph &g);


#endif //CPPPROJCT_ADEQUATE_SYNCHRONIZATION_H
