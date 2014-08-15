#!/bin/bash

if [ $# -ne 1 ] ; then
  echo "provide power_cap (e.g. 50)!" 
  exit 1
fi

cd output_$1     # goto power_cap output directory 

grep Reg * | grep tile > ../with_tile_$1.txt
#grep Reg * | grep parallel.c > ../wo_tile_$1.txt 

cd ..

#paste with_tile_$1.txt wo_tile_$1.txt  > with_wo_$1.txt 
cp with_tile_$1.txt with_wo_$1.txt


sort -k 4 -r -n with_wo_$1.txt > temp
mv temp with_wo_$1.txt 

cp with_wo_$1.txt data

gnuplot plot.gnpl
#gnuplot plot-png.gnpl
#gnuplot plot-time-energy-png.gnpl

mv power-time-compare.pdf $1.pdf
#mv power-time-compare.png $1-TE.png
gnuplot plot-time-energy.gnpl
mv power-time-compare.pdf $1-TE.pdf
