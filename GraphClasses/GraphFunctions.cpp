//
// Created by m4zz31 on 3/11/21.
//

#include "GraphFunctions.h"
#include <cassert>

#ifndef ASSERT
#define ASSERT
#endif

//void register_to_value(Graph &g){
//    auto vs = vertices(g);
//    int L = 0;
//#ifdef ASSERT
//    // Computing the value against which to compare
//    int tackled = 0;
//    for (auto v = vs.first; v != vs.second; ++v){
//        ++L;
//    }
//#endif
//#pragma omp parallel firstprivate(vs, L)
//{
//#ifdef ASSERT
//    int localtackled = 0;
//#endif
//    int N = omp_get_num_threads();
//    int i = omp_get_thread_num();
//    auto start = vs.first + (L/N) * i;
//    auto end = vs.first + (L/N) * (i + 1);
//    if (i == N-1) end += L % N;
//    for (auto v = start; v != end; ++v) {
//        printf("\nval and temp were: %f %f\n",g[*v].value, g[*v].temporal_register);
//        g[*v].value = g[*v].temporal_register;
//        printf("now are: %f %f\n\n",g[*v].value, g[*v].temporal_register);
//#ifdef ASSERT
//        localtackled++;
//#endif
//    }
//#ifdef ASSERT
//#pragma omp atomic update
//        tackled += localtackled;
//#endif
//}
//#ifdef ASSERT
//    assert(tackled == L);
//#endif
//};



void register_to_value(Graph &g){
    auto vs = vertices(g);
    int L = 0;

    // Computing the value against which to compare
    int tackled = 0;
    for (auto v = vs.first; v != vs.second; ++v){
        ++L;
    }
    int localtackled = 0;

    auto start = vs.first ;
    auto end = vs.second;
    for (auto v = start; v != end; ++v) {
        PRINTF_DBG("\nval and temp were: %f %f\n",g[*v].value, g[*v].temporal_register);
        g[*v].value = g[*v].temporal_register;
        PRINTF_DBG("now are: %f %f\n\n",g[*v].value, g[*v].temporal_register);
        localtackled++;
    }
#pragma omp atomic update
    tackled += localtackled;
};
