//
// Created by m4zz31 on 3/11/21.
//

#include "GraphFunctions.h"
#include <cassert>

#ifndef ASSERT
#define ASSERT
#endif
//

void update_neighbor_values(ReferenceContainer &REF){

    int N = (int) REF.NVtot;
    auto vs = vertices(*REF.p_g);
    int MYPROC = REF.p_ComHelper->WORLD_RANK[0];
#pragma omp parallel firstprivate(N, REF, MYPROC, vs)
{
    int Nthreads = (int) omp_get_num_threads();
    int myThread = (int) omp_get_thread_num();
    int begin = N/Nthreads * myThread;
    int end = N/Nthreads * (myThread + 1);
    if (myThread + 1 == Nthreads) end += N % Nthreads;
    // Update diagonal terms
    for (auto v = vs.first + begin; v != vs.first + end; ++v) {
        long ix = get(get(boost::vertex_index, *(REF.p_g)), *v);
        PRINTF_DBG("Updating central value for ix %ld processor %d: from %f to %f!\n",
                   ix, MYPROC, (*REF.p_IntHelper)[ix].centralValue, (*REF.p_g)[*v].value);
        (*REF.p_IntHelper)[ix].centralValue = (*REF.p_g)[*v].value;
    }
#pragma omp barrier
    // Propagate from diagonal terms
    for (int ix = begin; ix < end; ++ix){
        for (int j=0; j<(*REF.p_IntHelper)[ix].ixMap.size(); ++j){
            if (std::get<2>((*REF.p_IntHelper)[ix].ixMap[j]) == MYPROC){
                (*REF.p_IntHelper)[ix].neighborValues[j] = (*REF.p_IntHelper)[std::get<1>((*REF.p_IntHelper)[ix].ixMap[j])].centralValue;
            }
        }
    }
}
};


void register_to_value(Graph &g){
    auto vs = vertices(g);
    int L = 0;
#ifdef ASSERT
    // Computing the value against which to compare
    int tackled = 0;
    for (auto v = vs.first; v != vs.second; ++v){
        ++L;
    }
#endif
#pragma omp parallel firstprivate(vs, L)
{
#ifdef ASSERT
    int localtackled = 0;
#endif
    int N = omp_get_num_threads();
    int i = omp_get_thread_num();
    auto start = vs.first + (L/N) * i;
    auto end = vs.first + (L/N) * (i + 1);
    if (i == N-1) end += L % N;
    for (auto v = start; v != end; ++v) {
        PRINTF_DBG("\nval and temp were: %f %f. Vix: %lu proc: %d\n",
               g[*v].value,
               g[*v].temporal_register,
               get(get(boost::vertex_index, g), *v),
               boost::graph::distributed::process_id(g.process_group()));
        g[*v].value = g[*v].temporal_register;

#ifdef ASSERT
        localtackled++;
#endif
    }
#ifdef ASSERT
#pragma omp atomic update
        tackled += localtackled;
#endif
}
#ifdef ASSERT
    assert(tackled == L);
#endif
};


//
//void register_to_value(Graph &g){
//    auto vs = vertices(g);
//    int L = 0;
//
//    // Computing the value against which to compare
//    int tackled = 0;
//    for (auto v = vs.first; v != vs.second; ++v){
//        ++L;
//    }
//    int localtackled = 0;
//
//    auto start = vs.first ;
//    auto end = vs.second;
//    for (auto v = start; v != end; ++v) {
//        PRINTF_DBG("\nval and temp were: %f %f\n",g[*v].value, g[*v].temporal_register);
//        g[*v].value = g[*v].temporal_register;
//        PRINTF_DBG("now are: %f %f\n\n",g[*v].value, g[*v].temporal_register);
//        localtackled++;
//    }
//#pragma omp atomic update
//    tackled += localtackled;
//};
