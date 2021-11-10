//!nvcc -arch=sm_70 -o thread-and-block-idx 03-indices/01-thread-and-block-idx.cu -run


// set NThreads as a multiple of 32 for performance optimization of CPU instructions

// max threads per block are  1024

//cudaFree and   cudaMallocManaged(&a, size) instead of malloc and free


//make a directory of kernels


// profile with !nsys profile --stats=true ./iteratively-optimized-vector-add

// assign #StreamProcessors to #Blocks for perofrmance opt. (load balance)
