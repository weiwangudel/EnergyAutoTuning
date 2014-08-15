#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

rm NoBarrier_TEP_avg
rm Barrier_TEP_avg
rm Original_TEP_avg
 
for i in `seq 40 5 100` 0
do
   grep Parallel  NoBarrier_$i  |sort -k 5 -r -n | tail -n 5 | cut -d" " -f5,9,13 >> NoBarrier_TEP_4x
   grep Parallel  NoBarrier_$i  |sort -k 5 -r -n | tail -n 5 | cut -d" " -f5,9,13 > temp
   
   awk '{t+=$1;e+=$2;p+=$3} END{print '$i' " " t/NR " " e/NR " " p/NR}' temp >> NoBarrier_TEP_avg
   rm temp
   
   grep Parallel  Barrier_$i  |sort -k 5 -r -n | tail -n 5 | cut -d" " -f5,9,13 >> Barrier_TEP_4x
   grep Parallel  Barrier_$i  |sort -k 5 -r -n | tail -n 5 | cut -d" " -f5,9,13 > temp
   
   awk '{t+=$1;e+=$2;p+=$3} END{print '$i' " " t/NR " " e/NR " " p/NR}' temp >> Barrier_TEP_avg
   rm temp

   grep Parallel  Original_$i  |sort -k 5 -r -n | tail -n 5 | cut -d" " -f5,9,13 >> Original_TEP_4x
   grep Parallel  Original_$i  |sort -k 5 -r -n | tail -n 5 | cut -d" " -f5,9,13 > temp
   
   awk '{t+=$1;e+=$2;p+=$3} END{print '$i' " " t/NR " " e/NR " " p/NR}' temp >> Original_TEP_avg
   rm temp

done 
