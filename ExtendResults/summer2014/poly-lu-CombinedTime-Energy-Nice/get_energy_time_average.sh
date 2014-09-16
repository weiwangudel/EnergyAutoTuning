#!/bin/bash

for i in `seq 6 1 16`; do grep Reg lu-dc-$i.output  |sort -k 4 -r -n    |awk '{t+=$4;e+=$8;p+=$12}END{print '$i'/16*100"%\t" t/NR "\t" e/NR "\t" p/NR } ' ; done
