#!/bin/bash

if [ $# -ne 1 ] ; then
  echo "provide loop #  as first parameter"
  exit 1
fi 

for dc in `seq 2700000 -100000 2100000` `seq 1900000 -100000 1200000`  
do
     grep "Loop $1" 16-$dc.output  |grep Memory > temp
   
     awk '{t+=$6;e+=$10;p+=$14} END{print '$dc' " " t/NR " " e/NR " " p/NR}' temp
          

     rm temp
done
