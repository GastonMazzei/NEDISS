//
// Created by m4zz31 on 6/11/21.
//

#include "print_init.h"
#include "adequate_synchronization.h"

void print_init(int rank){
    if (rank == 0) {
        std::cout << init_msg << std::endl;
    }
};