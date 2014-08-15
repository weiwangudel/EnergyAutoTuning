#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

rm Original_TEP_4x
rm WithDutyCycle_TEP_4x
rm orig_TEP_avg
rm withDutyCycle_TEP_avg
 
for i in `seq 40 5 100` 0
do
   grep socket  Original_$i  |tail -n 4 | cut -d" " -f4,8,12 >> Original_TEP_4x
   grep socket  Original_$i  |tail -n 4 | cut -d" " -f4,8,12 > temp
   
   awk '{t+=$1;e+=$2;p+=$3} END{print '$i' " " t/NR " " e/NR " " p/NR}' temp >> orig_TEP_avg
   rm temp
   
   grep socket  Barrier_$i  |tail -n 4 | cut -d" " -f4,8,12 >> WithBarrier_TEP_4x
   grep socket  Barrier_$i  |tail -n 4 | cut -d" " -f4,8,12 > temp
   
   awk '{t+=$1;e+=$2;p+=$3} END{print '$i' " " t/NR " " e/NR " " p/NR}' temp >> withBarrier_TEP_avg
   rm temp
done 
