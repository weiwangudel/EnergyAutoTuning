#!/bin/bash

#if [ $# -ne 1 ] ; then
#  echo "provide 1 (time), 2 (energy), 3 (power) as first parameter"
#  exit 1
#fi 

#for freq in 2700000 2400000 2100000 1800000 1500000 1200000

for i in 4 5 6 7 8 9
do

   sed 's/Wei/'$i'/g' draw-template.gnpl  > draw.gnpl
   gnuplot draw.gnpl
   mv output.pdf $i.pdf
done
