#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

#for freq in 2701000 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000
for freq in 2700000 2400000 2100000 1800000 1500000 
do
#   for num in 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
  # for num in 16 15 14 13 12 
#   for num in 1 2 4 8 16 
   for num in 16 15 14 12 10 8

   do 
#     echo "freq$freq-num$num"
     for variants in `ls *freq$freq-num$num.XLoutput`
     do
       grep socket $variants  |cut -d " " -f4,8,12,18,23 > temp
       echo $variants
       #awk '{t+=$1;e+=$2;p+=$3} END{print "data" " " '$freq' " " '$num' " "  t/NR " " e/NR " " p/NR}' temp
       awk '{t+=$1;e+=$2;p+=$3;t1+=$4;t2+=$5} END{print "data" " " t/NR " " e/NR " " p/NR " " t1/NR " " t2/NR}' temp
       rm temp
     done
   done
   echo ""
done
   
