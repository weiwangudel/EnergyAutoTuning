#!/bin/bash

if [ $# -ne 1 ] ; then
  echo "provide hybrid output file as first parameter"
  exit 1
fi 

#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
       #grep Parallel 
       grep Parallel   $1 |awk '{print $8}' > parallel
       grep App        $1|grep socket |awk '{print $8}' > App
       grep Reg        $1|awk '{print $8}' > IO
       grep "Loop 13 " $1|grep Power |awk '{print $10}' > 13
       grep "Loop 14 " $1|grep Power |awk '{print $10}' > 14
       grep "Loop 16 " $1|grep Power |awk '{print $10}' > 16

       paste App IO parallel 13 14 16 > pA345
       rm App IO parallel 13 14 16
       #sort -k 2 -r -n pA345-$dc | awk '{a+=$1;io+=$2;p+=$3;c+=$4;d+=$5;e+=$6} END {print a/NR " " io/NR " " p/NR " " c/NR " " d/NR " " e/NR}' 
       cat pA345 | awk '{a+=$1;io+=$2;p+=$3;c+=$4;d+=$5;e+=$6} END {print a/NR "\t" io/NR "\t" p/NR "\t" c/NR "\t" d/NR "\t" e/NR}' 
       rm pA345
   
	
