#!/bin/bash

cat only_top100_smartfuse_tiling_size_for_maxfuse_as_well.csv | while read i j k 
do
   echo "$i $j $k" > tile.sizes
   cp Selected-rose_lulesh-body-$i"x"$j"x"$k.c.exe.output    Selected-100/
done 
