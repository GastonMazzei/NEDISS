#include <iostream>
#include <vector>
#include <iostream>
#include <random>
#include <omp.h>
#include <algorithm>
#include <iterator>
#include <chrono>

void time_L_executions(int L){
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int j=0; j<100000; j++) {
        for (int i = 0; i < L; i++) {
            //std::cout << "Hello, World!" << std::endl;
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count()/100000 << std::endl;
}

void parallel_time_L_executions(int L){
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int j=0; j<100000; j++) {
# pragma omp parallel for
        for (int i = 0; i < L; i++) {
            //std::cout << "Hello, World!" << std::endl;
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count()/100000 << std::endl;
}

// TO JOIN TIMERS ITERATING OVER THEM:
//
//for (int L=0.5e2; L < 2.5e3; L = L + 1e2) {
//std::cout << L << std::endl;
//time_L_executions(L);
//parallel_time_L_executions(L);
//}