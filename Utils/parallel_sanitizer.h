//
// Created by m4zz31 on 6/11/21.
//

#ifndef CPPPROJCT_PARALLEL_SANITIZER_H
#define CPPPROJCT_PARALLEL_SANITIZER_H

#include "global_standard_messages.h"
#include "error.h"

void check_nested_status(){
    char const* tmp = getenv( "OMP_NESTED" );
    if (tmp==NULL){
        error_report(nested_threads);
    }
    std::string s( tmp );
    if ((s!="true") && (s!="TRUE")){
        error_report(nested_threads);
    }
};



#endif //CPPPROJCT_PARALLEL_SANITIZER_H
