#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
for dc in `seq 6 1 16` 
do
       grep Reg  covariance-dc-$dc.output  |awk '{print $4}'  > parallel

       sort -k 2 -r -n parallel | awk '{e+=$1;} END {print '$dc' "\t" '$dc'/16*100 "%\t" e/NR }' 
       rm parallel
done
   
