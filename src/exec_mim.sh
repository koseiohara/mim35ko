#!/bin/bash

#PBS -q tqueue
#PBS -N MIM_ALL
#PBS -j oe
#PBS -l nodes=1:ppn=1

cd /mnt/hail8/kosei/mim/JRA3Q/mymim35ALL/src/
ulimit -s unlimited

NOW=$(date "+%Y%m%d_%H%M%S")
INI=1970
FIN=2022
RESULT_FILE="../output/result_${INI}_${FIN}_${NOW}.txt"
./MIM_ALL3Q >& ${RESULT_FILE}

