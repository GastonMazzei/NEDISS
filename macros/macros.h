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

#if !defined(VERTEXVAL_REQUEST_FLAG)
#define VERTEXVAL_REQUEST_FLAG 25001
#endif

#if !defined(EDGEVAL_REQUEST_FLAG)
#define EDGEVAL_REQUEST_FLAG 25000
#endif

#if !defined(OFFSET)
#define OFFSET 12500
#endif


#endif //CPPPROJCT_MACROS_H
