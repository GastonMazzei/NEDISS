#!/bin/bash
START=$(date +%s)
mpirun -oversubscribe -x NNODES=50 -x TEST=1 -x SEED=21234 -x OMP_THREAD_LIMIT=6  -x  OMP_NESTED=true -np 3  cmake-build-debug/cppprojct > tmp/current.log
rm  tmp/current.log
# your logic ends here
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
