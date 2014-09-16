#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
for freq in 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000

do
       grep App       freq-$freq.output |grep socket |awk '{print $8}' > App
       grep Parallel  freq-$freq.output  |awk '{print $8}' > parallel
       grep "Loop 2 " freq-$freq.output |grep Power |awk '{print $10}' > 2
       grep "Loop 7 " freq-$freq.output |grep Power |awk '{print $10}' > 7
       grep "Loop 0 " freq-$freq.output |grep Power |awk '{print $10}' > 0 
       grep "Loop 1 " freq-$freq.output |grep Power |awk '{print $10}' > 1
       grep "Loop 3 " freq-$freq.output |grep Power |awk '{print $10}' > 3
       grep "Loop 8 " freq-$freq.output |grep Power |awk '{print $10}' > 8

       paste App parallel 2 7 0 1 3 8 > Ap270138-$freq
       sort -k 2 -r -n Ap270138-$freq | awk '{a+=$1;p+=$2;b+=$3;c+=$4;d+=$5;e+=$6;f+=$7;g+=$8;} END {print "'$freq'" " "  a/NR " " p/NR " " b/NR " " c/NR " " d/NR " " e/NR " " f/NR " "g/NR }' 
       rm Ap270138-$freq
       rm 2 7 0 1 3 8 App parallel
done
   
