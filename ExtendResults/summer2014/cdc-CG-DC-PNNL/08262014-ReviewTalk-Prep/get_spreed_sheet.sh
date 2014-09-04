#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

#for freq in 2701000 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000
for freq in 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000
#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
do
   for num in 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
#   for num in 16 15 14 12 10 8 
#   for num in 1 2 4 8 16 

   do 
  #   grep  line-99 freq-$freq-num$num.output  |cut -d " " -f6,10,14 > temp
     grep Parallel freq-$freq-num$num.output  |cut -d " " -f5,9,13 > temp
   
     awk '{t+=$1;e+=$2;p+=$3} END{print '$freq' " " '$num' " "  t/NR " " e/NR " " p/NR}' temp
     rm temp
   done
   echo ""
done
   
