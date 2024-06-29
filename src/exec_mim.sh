#!/bin/bash

#PBS -q tqueue
#PBS -N MIM_ALL_FAST
#PBS -j oe
#PBS -l nodes=1:ppn=1

START=$(date "+%s")

cd /mnt/jet11/kosei/mim/mim35ko/src/
ulimit -s unlimited

NOW=$(date "+%Y%m%d_%H%M%S")
INI=2000
FIN=2000
RESULT="../output/result_${INI}_${FIN}_${NOW}.txt"

NAMELIST="../nml/input_JRA55_2000_2000.nml"

./MIM_ALL3Q < ${NAMELIST} >& ${RESULT}

END=$(date "+%s")


DIFF_SEC=$(expr ${END} - ${START})
DIFF_MIN=$(expr ${DIFF_SEC} / 60)
DIFF_SEC=$(expr ${DIFF_SEC} - ${DIFF_MIN} \* 60)
DIFF_HR=$(expr ${DIFF_MIN} / 60)
DIFF_MIN=$(expr ${DIFF_MIN} - ${DIFF_HR} \* 60)
echo " " >> ${RESULT}
echo "ELAPS : ${DIFF_HR}hr ${DIFF_MIN}min ${DIFF_SEC}sec" >> ${RESULT}

cat ${NAMELIST} >> ${RESULT}

