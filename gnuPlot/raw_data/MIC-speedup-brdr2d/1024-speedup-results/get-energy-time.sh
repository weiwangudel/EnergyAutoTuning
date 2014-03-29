#/bin/bash

grep "Region" base_power* 1024-102* |awk '(NR == 1) {v1=$4; v2=$8;} (NR>1) {print v1/$4, 1-$8/v2}' 
