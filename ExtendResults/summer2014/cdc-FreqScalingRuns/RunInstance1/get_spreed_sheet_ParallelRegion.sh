#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

for freq in 2700000 2400000 2100000 1800000 1500000 1200000
do
   for num in 1 2 4 8 16 
   do 
     grep Parallel freq-$freq-num$num.output  |cut -d " " -f5,9,13 > temp
   
     awk '{t+=$1;e+=$2;p+=$3} END{print '$freq' " " '$num' " "  t/NR " " e/NR " " p/NR}' temp
     rm temp
   done
   echo ""
done
   
