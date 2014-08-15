#!/bin/bash

if [ $# -ne 2 ] ; then
  echo "provide power_cap (e.g. 50 benchmark)!" 
  exit 1
fi

cd output_$1     # goto power_cap output directory 

grep Reg *  > ../with_powercap_$1.txt

cd ..
cd output_1000
grep Reg * > ../with_powercap_1000.txt

cd ..

paste with_powercap_1000.txt  with_powercap_$1.txt > with_wo_$1.txt 
#
sort -k 4 -r -n with_wo_$1.txt > temp
mv temp with_wo_$1.txt 

cp with_wo_$1.txt data

gnuplot plot.gnpl
gnuplot plot-png.gnpl
#
mv power-time-compare.pdf $2-$1.pdf
mv power-time-compare.png $2-$1.png
