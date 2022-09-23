#!/bin/bash
range=300


for i in `seq 0 $range`;
do
    echo $i     # iteration print
    gcc rtn_simple.c -o rtn -lm -std=c99 -lfftw3l && ./rtn -t $(cat t0_out.txt)     # compilation and execution of the rtn code with the new tc[0] parameter extracted from a txt file
done