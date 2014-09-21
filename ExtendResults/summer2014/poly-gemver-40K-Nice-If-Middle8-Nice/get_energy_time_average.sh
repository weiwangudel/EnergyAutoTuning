#!/bin/bash

if [ $# -lt 1 ] ; then
  echo " $0 polybench-name head tail "
  exit  1
fi

for i in `seq 6 1 16`; do grep Reg $1-dc-$i.output  |sort -k 4 -r -n  |head -n $2 |tail -n $3  |awk '{t+=$4;e+=$8;p+=$12}END{print '$i'/16*100"%\t" t/NR "\t" e/NR "\t" p/NR } ' ; done
