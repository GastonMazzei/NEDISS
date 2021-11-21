#!/bin/bash
START=$(date +%s)
./run.sh > tmp/current.log
rm  tmp/current.log
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF second(s)"
