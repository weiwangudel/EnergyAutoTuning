#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
for dc in `seq 1 1 16` 
do
   for num in 16 
   do 
       #grep Parallel 
       grep Parallel   16threads-$dc.output  |awk '{print $4}' > parallel
       grep App        16threads-$dc.output |grep socket |awk '{print $4}' > App
       grep Reg        16threads-$dc.output |awk '{print $4}' > IO
       grep "Loop 13 " 16threads-$dc.output |grep Power |awk '{print $6}' > 13
       grep "Loop 14 " 16threads-$dc.output |grep Power |awk '{print $6}' > 14
       grep "Loop 16 " 16threads-$dc.output |grep Power |awk '{print $6}' > 16

       paste App IO parallel 13 14 16 > pA345-$dc
       sort -k 2 -r -n pA345-$dc | awk '{a+=$1;io+=$2;p+=$3;c+=$4;d+=$5;e+=$6} END {print '$dc'/16*100 "% " a/NR " " io/NR " " p/NR " " c/NR " " d/NR " " e/NR}' 
   done
done
   
