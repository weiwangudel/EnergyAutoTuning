#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

for dc in `seq 16 -1 1` 
do
     grep Parallel 16threads-$dc.output  |cut -d " " -f5,9,13 > temp
   
     awk '{t+=$1;e+=$2;p+=$3} END{print '$dc' " " '$num' " "  t/NR " " e/NR " " p/NR}' temp
          

     rm temp
done
