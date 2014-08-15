#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Which Loop?"
  exit  1
fi

    awk '{print $1 " " $3 " " $2*$3 }'  data  > temp
    
    minE=$(sort -k 2 -r -n temp |tail -n 1 | cut -d" " -f2)
    minED=$(sort -k 3 -r -n temp |tail -n 1 | cut -d" " -f3)
    
    
   awk '{print $1/1000000 " "  $2/'$minE' " " $3/'$minED' }' temp > gnuplot.data 
   tac gnuplot.data > temp
   mv temp gnuplot.data
   gnuplot draw.gnpl
   mv output.pdf $1.pdf
