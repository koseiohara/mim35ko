#!/bin/bash

#PBS -q tqueue
#PBS -N MIM_ALL3Q_FAST
#PBS -j oe
#PBS -l nodes=1:ppn=1

START=$(date "+%s")

cd /mnt/jet11/kosei/mim/mim35ko/src/
ulimit -s unlimited

NOW=$(date "+%Y%m%d_%H%M%S")
INI=1975
FIN=2022
RESULT_FILE="../output/result_${INI}_${FIN}_${NOW}.txt"
./MIM_ALL3Q >& ${RESULT_FILE}

END=$(date "+%s")


DIFF_SEC=$(expr ${END} - ${START})
DIFF_MIN=$(expr ${DIFF_SEC} / 60)
DIFF_SEC=$(expr ${DIFF_SEC} - ${DIFF_MIN} \* 60)
DIFF_HR=$(expr ${DIFF_MIN} / 60)
DIFF_MIN=$(expr ${DIFF_MIN} - ${DIFF_HR} \* 60)
echo " " >> ../test.txt
echo "ELAPS : ${DIFF_HR}hr ${DIFF_MIN}min ${DIFF_SEC}sec" >> ${RESULT_FILE}

