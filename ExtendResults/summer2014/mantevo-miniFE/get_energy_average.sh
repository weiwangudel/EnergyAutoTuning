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
       grep App       dc-300-$dc.output |grep socket |awk '{print $8}'  > App
       grep Parallel  dc-300-$dc.output  |awk '{print $8}'  > parallel
       grep "Loop 0 " dc-300-$dc.output |grep Power |awk '{print $10}' > 0
       grep "Loop 1 " dc-300-$dc.output |grep Power |awk '{print $10}' > 1 
       grep "Loop 2 " dc-300-$dc.output |grep Power |awk '{print $10}' > 2 
       grep "Loop 3 " dc-300-$dc.output |grep Power |awk '{print $10}' > 3 

       paste App parallel 0 1 2 3 > Ap0123-$dc
       sort -k 2 -r -n Ap0123-$dc | awk '{a+=$1;p+=$2;b+=$3;c+=$4;d+=$5;e+=$6;} END {print '$dc' "\t" '$dc'/16*100 "%\t" a/NR "	" p/NR "	" b/NR "	" c/NR "	" d/NR "	" e/NR }' 
       rm Ap0123-$dc
       rm 0 1 2 3 App parallel
   done
done
   
