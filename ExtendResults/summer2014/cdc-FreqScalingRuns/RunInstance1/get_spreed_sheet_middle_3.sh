#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

for freq in 2700000 2400000 2100000 1800000 1500000 1200000
do
   for num in 1 2 4 8 16 
   do 
     grep socket freq-$freq-num$num.output  |sort -k 4 -r -n | tail -n 3 | cut -d " " -f4,8,12 > temp
   
     awk '{t+=$1;e+=$2;p+=$3} END{print '$freq' " " '$num' " "  t/NR " " e/NR " " p/NR}' temp
     rm temp
   done
   echo ""
done
   
