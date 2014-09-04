#!/bin/bash

if [ $# -ne 1 ] ; then
  echo "provide hybrid output file as first parameter"
  exit 1
fi 

#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
       #grep Parallel 
       grep App       $1 |grep socket |awk '{print $8}' > App
       grep Parallel  $1  |awk '{print $8}' > parallel
       grep "Loop 2 " $1|grep Power |awk '{print $10}' > 2
       grep "Loop 7 " $1|grep Power |awk '{print $10}' > 7
       grep "Loop 0 " $1|grep Power |awk '{print $10}' > 0 
       grep "Loop 1 " $1|grep Power |awk '{print $10}' > 1
       grep "Loop 3 " $1|grep Power |awk '{print $10}' > 3
       grep "Loop 8 " $1|grep Power |awk '{print $10}' > 8

       paste App parallel 2 7 0 1 3 8 > Ap270138
       sort -k 2 -r -n Ap270138 |head -n 5 | awk '{a+=$1;p+=$2;b+=$3;c+=$4;d+=$5;e+=$6;f+=$7;g+=$8;} END {print a/NR "\t" p/NR "\t" b/NR "\t" c/NR "\t" d/NR "\t" e/NR "\t" f/NR "\t"g/NR }' 
       cat Ap270138
       rm Ap270138
       rm 2 7 0 1 3 8 App parallel


   
	
