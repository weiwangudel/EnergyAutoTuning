#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
for freq in 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000
do
   for num in 16 
   do 
       #grep Parallel 
       grep Parallel freq-$freq-num16.output  |awk '{print $8}' > parallel
       grep App freq-$freq-num16.output |grep socket |awk '{print $8}' > App
       grep Reg freq-$freq-num16.output |awk '{print $8}' > IO
       grep "Loop 13 " freq-$freq-num16.output |grep Power |awk '{print $10}' > 13
       grep "Loop 14 " freq-$freq-num16.output |grep Power |awk '{print $10}' > 14
       grep "Loop 16 " freq-$freq-num16.output |grep Power |awk '{print $10}' > 16

       paste App IO parallel 13 14 16 > pA345-$dc
       sort -k 2 -r -n pA345-$dc | awk '{a+=$1;io+=$2;p+=$3;c+=$4;d+=$5;e+=$6} END {print a/NR " " io/NR " " p/NR " " c/NR " " d/NR " " e/NR}' 

   done
done
   
