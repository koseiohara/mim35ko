#!/bin/bash

#PBS -q tqueue
#PBS -N MIM_FAST_552000
#PBS -j oe
#PBS -l nodes=1:ppn=1

START=$(date "+%s")

cd /mnt/jet11/kosei/mim/mim35ko/src/
ulimit -s unlimited

NOW=$(date "+%Y%m%d_%H%M%S")
DATA="JRA55"
INI=2000
FIN=2000
RESULT="../output/${DATA}/result_${INI}_${FIN}_${NOW}.txt"

#NAMELIST="../nml/input_JRA55_2000_2000.nml"
NAMELIST="../nml/input_${DATA}_${INI}_${FIN}_fast.nml"

./MIM < ${NAMELIST} >& ${RESULT}

END=$(date "+%s")


DIFF_SEC=$(expr ${END} - ${START})
DIFF_MIN=$(expr ${DIFF_SEC} / 60)
DIFF_SEC=$(expr ${DIFF_SEC} - ${DIFF_MIN} \* 60)
DIFF_HR=$(expr ${DIFF_MIN} / 60)
DIFF_MIN=$(expr ${DIFF_MIN} - ${DIFF_HR} \* 60)
echo " " >> ${RESULT}
echo "ELAPS : ${DIFF_HR}hr ${DIFF_MIN}min ${DIFF_SEC}sec" >> ${RESULT}

echo -e "\n\n" >> ${RESULT}
cat ${NAMELIST} >> ${RESULT}

