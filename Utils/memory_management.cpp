//
// Created by m4zz31 on 1/11/21.
//

#include "memory_management.h"

void clear_vectors(std::vector<double> a,
                std::vector<double> b,
                 std::vector<double> c){
    a.clear();
    a.shrink_to_fit();
    b.clear();
    b.shrink_to_fit();
    c.clear();
    c.shrink_to_fit();
};