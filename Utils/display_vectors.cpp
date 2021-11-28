//
// Created by m4zz31 on 28/11/21.
//
#include "display_vectors.h"

void display(std::vector<double> a){
    for (int i=0; i<a.size();++i){
        std::cout << a[i] << " , ";
    }
    std::cout << std::endl;
}
