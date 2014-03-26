#!/bin/bash

cat only_top100_smartfuse_tiling_size_for_maxfuse_as_well.csv | while read i j k 
do
   echo "$i $j $k" > tile.sizes
    timeout 1200 polyopt -lm -std=c99 lulesh-body.c --polyopt-safe-math-func --polyopt-pluto-scalpriv  --polyopt-pluto-fuse-maxfuse --polyopt-pluto-parallel --polyopt-pluto-tile --polyopt-pluto-prevector --polyopt-generate-pragmas  
    if [ -f "rose_lulesh-body.c" ] 
    then
        cp rose_lulesh-body.c maxfuse_rose_lulesh-body-$i"x"$j"x"$k.c 
        rm rose_lulesh-body.c
    else
        echo "Polyopt was not successful in transform with tile size $i x $j x $k"
    fi

done 
