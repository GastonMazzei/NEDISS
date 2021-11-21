//
// Created by m4zz31 on 10/11/21.
//

#ifndef CPPPROJCT_MACROS_H
#define CPPPROJCT_MACROS_H

// lEGACY definition of VERBOSE = true or false here... should migrate!
#if defined(VERBOSE)
#define PRINTF_DBG printf
#else
#define PRINTF_DBG(...)
#define VERBOSE false
#endif

#if !defined(INTEGRATION_PLACEHOLDER)
#define INTEGRATION_PLACEHOLDER 1.010101
#endif

#if !defined(VERTEXVAL_REQUEST_FLAG)
#define VERTEXVAL_REQUEST_FLAG 25001
#endif

#if !defined(EDGEVAL_REQUEST_FLAG)
#define EDGEVAL_REQUEST_FLAG 25000
#endif

#if !defined(K1_REQUEST)
#define K1_REQUEST 24999
#endif

#if !defined(K1_ANSWER)
#define K1_ANSWER 24998
#endif

#if !defined(K2_REQUEST)
#define K2_REQUEST 24997
#endif

#if !defined(K2_ANSWER)
#define K2_ANSWER 24996
#endif

#if !defined(K3_REQUEST)
#define K3_REQUEST 24995
#endif

#if !defined(K3_ANSWER)
#define K3_ANSWER 24994
#endif

#if !defined(K4_REQUEST)
#define K4_REQUEST 24993
#endif

#if !defined(K4_ANSWER)
#define K4_ANSWER 24992
#endif



#if !defined(OFFSET)
#define OFFSET 12500
#endif

//#define ASSERTS

#endif //CPPPROJCT_MACROS_H
