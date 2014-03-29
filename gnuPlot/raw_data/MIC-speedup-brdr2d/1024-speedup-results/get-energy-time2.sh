#/bin/bash

grep "Region" base_power* 1024-512* |awk '(NR == 1) {v1=$4; v2=$8;} (NR>1) {print v1/$4, v2/$8}' 
