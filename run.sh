

mpiexec -x OMP_THREAD_LIMIT=4  -x  OMP_NESTED=true -n 2  cmake-build-debug/cppprojct
#mpiexec -x OMP_THREAD_LIMIT=$1  -x  OMP_NESTED=true -n $2  cmake-build-debug/cppprojct
