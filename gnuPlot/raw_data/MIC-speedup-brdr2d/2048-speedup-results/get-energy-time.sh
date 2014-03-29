#/bin/bash

grep "Region" *.output |awk '(NR == 1) {v1=$4; v2=$8;} (NR>1) {print v1/$4, 1-$8/v2}' 
