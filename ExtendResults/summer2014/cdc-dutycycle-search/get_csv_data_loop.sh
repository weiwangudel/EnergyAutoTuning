#!/bin/bash

if [ $# -ne 1 ] ; then
  echo "provide loop #  as first parameter"
  exit 1
fi 

for dc in `seq 16 -1 1`
do
     grep "Loop $1" 16-$dc.output  |grep Memory > temp
   
     awk '{t+=$6;e+=$10;p+=$14} END{print '$dc' " " t/NR " " e/NR " " p/NR}' temp
          

     rm temp
done
