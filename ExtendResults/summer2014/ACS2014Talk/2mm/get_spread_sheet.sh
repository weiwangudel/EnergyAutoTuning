#!/bin/bash


for cap in 0 40
do
  rm $cap
  cd output_$cap
   
  for i in `cat ../2mm-names` 
  do
      echo $i > name 
      grep Reg $i* |cut -d" " -f4,8,12 > temp
      awk '{t+=$1;e+=$2;p+=$3} END {print t/NR " " e/NR " " p/NR " " }' temp > temp2
      
      paste -d" " name temp2 >> ../$cap
      rm name
      rm temp 
      rm temp2
  done 
  cd ..
done 


paste -d" " 0  40 > with_wo_40.txt 

sort -k 2 -r -n with_wo_40.txt > temp
mv temp with_wo_40.txt 

cp with_wo_40.txt data

gnuplot plot.gnpl
#gnuplot plot-png.gnpl
#
mv power-time-compare.pdf 2mm-40.pdf
#mv power-time-compare.png 2mm-40.png
