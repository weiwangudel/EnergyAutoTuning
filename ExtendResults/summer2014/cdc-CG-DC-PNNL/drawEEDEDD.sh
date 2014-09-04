#!/bin/bash

    awk '{print $1 " " $3 " " $2*$3 }'  data  > temp
    
    minE=$(sort -k 2 -r -n temp |tail -n 1 | cut -d" " -f2)
    minED=$(sort -k 3 -r -n temp |tail -n 1 | cut -d" " -f3)
    
    
   awk '{print $1 " "  $2/'$minE' " " $3/'$minED' }' temp > gnuplot.data 
   tac gnuplot.data > temp
   mv temp gnuplot.data
   gnuplot draw.gnpl
   mv output.pdf dutycycle-parallel.pdf
