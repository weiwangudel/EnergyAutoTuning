#!/bin/bash

if [ $# -lt 1 ] ; then 
  echo "$0 polybench-name"
  exit 1
fi

for dc in `seq 16 -1 12`
do 
  for num in `seq 16 -2 2`
  #`seq 2 2 16` 
  do 
      grep "Loop 3 " $1-dc-$dc-num-$num.output |sort -k 6 -r -n |head -n 7 |tail -n 5 |awk '{t+=$6;e+=$10;p+=$14;} END {print '$dc'/16*100 "% " '$num' " " t/NR  " "e/NR " " p/NR}'  
  done
  echo " "
done
