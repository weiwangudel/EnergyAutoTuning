#!/bin/bash

if [ $# -lt 1 ] ; then 
  echo "$0 polybench-name"
  exit 1
fi

for dc in `seq 16 -1 12`
do 
  for num in `seq 2 2 16` 
  do 
      grep App $1-icc-dc-$dc-num-$num.output |sort -k 4 -r -n |head -n 7 |tail -n 5 |awk '{t+=$4;e+=$8;p+=$12;} END {print '$dc'/16*100 "% " '$num' " " t/NR  " "e/NR " " p/NR}'  
  done
done
