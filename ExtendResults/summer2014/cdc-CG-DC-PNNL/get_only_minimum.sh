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
       grep Parallel   16threads-$dc.output  |awk '{print $8}' > parallel
       grep App        16threads-$dc.output |grep socket |awk '{print $8}' > App
       grep "Loop 13 " 16threads-$dc.output |grep Power |awk '{print $10}' > 13
       grep "Loop 14 " 16threads-$dc.output |grep Power |awk '{print $10}' > 14
       grep "Loop 16 " 16threads-$dc.output |grep Power |awk '{print $10}' > 16

       paste App parallel 13 14 16 > pA345-$dc
       sort -k 2 -r -n pA345-$dc | tail -n 8 |head -n 6  |awk '{a+=$1;p+=$2;c+=$3;d+=$4;e+=$5} END {print a/NR " " p/NR " " c/NR " " d/NR " " e/NR}' 
   done
   echo ""
done
   
