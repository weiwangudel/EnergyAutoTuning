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
       grep App       dc-120-$dc.output |grep socket |awk '{print $4}' > App
       grep Parallel  dc-120-$dc.output  |awk '{print $4}' > parallel
       grep "Loop 2 " dc-120-$dc.output |grep Power |awk '{print $6}' > 2
       grep "Loop 7 " dc-120-$dc.output |grep Power |awk '{print $6}' > 7
       grep "Loop 0 " dc-120-$dc.output |grep Power |awk '{print $6}' > 0 
       grep "Loop 1 " dc-120-$dc.output |grep Power |awk '{print $6}' > 1
       grep "Loop 3 " dc-120-$dc.output |grep Power |awk '{print $6}' > 3
       grep "Loop 8 " dc-120-$dc.output |grep Power |awk '{print $6}' > 8

       paste App parallel 2 7 0 1 3 8 > Ap270138-$dc
       sort -k 2 -r -n Ap270138-$dc | awk '{a+=$1;p+=$2;b+=$3;c+=$4;d+=$5;e+=$6;f+=$7;g+=$8;} END {print "'$dc'" " "  a/NR " " p/NR " " b/NR " " c/NR " " d/NR " " e/NR " " f/NR " "g/NR }' 
       rm Ap270138-$dc
       rm 2 7 0 1 3 8 App parallel
   done
done
   
