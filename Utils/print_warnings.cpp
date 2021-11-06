//
// Created by m4zz31 on 6/11/21.
//

#include "print_warnings.h"
#include "adequate_synchronization.h"


void print_warnings(int rank){
    if (rank == 0) {
        std::cout << init_warning << std::endl;
    }
};