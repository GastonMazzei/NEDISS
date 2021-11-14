//
// Created by m4zz31 on 13/11/21.
//
#include "certify_compliant.h"
#include "global_standard_messages.h"
#include <iostream>
#include "error.h"

void certify_tagMax_compliant(int TAGMAX_FLAG,
                              int MAXTAG,
                              const int &NNODES,
                              int NPROCS,
                              int MAX_USED_TAG){
    if (TAGMAX_FLAG == 1){
        if ((MAX_USED_TAG > MAXTAG) && (NNODES/NPROCS + 1 > MAXTAG)){
            std::string local_message = architecture_and_maxtag +
                                       "(1) the current MPI implementation does not support so many nodes per processor because of the maximum TAG value and how this program is designed." +
                                       " The only solution is raising the number of Processors to be >=" +
                                       std::to_string(NNODES/MAXTAG + 1) +
                                       " or lowering the number of nodes below " +
                                       std::to_string(NPROCS*MAXTAG - 1) + ". " +
                                       "(2) the compilation-defined constant 'VERTEXVAL_REQUEST_FLAG' currently exceeds the maximum allowed tag value. Please raise it to <=" +
                                       std::to_string(MAXTAG) +
                                       ", and if you chose to make it not equal but lower recompute the previous upper bounds (i.e. 1) for this tag value.";
            error_report(local_message);
        } else if (NNODES/NPROCS > MAXTAG){
            std::string local_message = architecture_and_maxtag +
                                       " the current MPI implementation does not support so many nodes per processor because of the maximum TAG value and how this program is designed." +
                                       " The only solution is raising the number of Processors to be >=" +
                                       std::to_string(NNODES/MAXTAG + 1) +
                                       " or lowering the number of nodes below " +
                                       std::to_string(NPROCS*MAXTAG - 1) + ". " ;
            error_report(local_message);
        } else if (MAX_USED_TAG > MAXTAG){
            std::string local_message = architecture_and_maxtag +
                                       "the compilation-defined constant 'VERTEXVAL_REQUEST_FLAG' currently exceeds the maximum allowed tag value. Please raise it to <=" +
                                       std::to_string(MAXTAG);
            error_report(local_message);
        }
    } else {
        std::string local_message =  "[warning][severe] The program cannot determine if there is a critical error as MPI implementation did not provide with maximum tag value";
        std::cout << local_message << std::endl;
        std::cout << std::flush;
    }
}
