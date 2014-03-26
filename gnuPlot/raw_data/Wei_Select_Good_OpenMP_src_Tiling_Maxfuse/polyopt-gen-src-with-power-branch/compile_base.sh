#!/bin/bash

for i in `ls lulesh*omp*.c `
do
  sed -i 's/1.0e-6/1.0e-4/g' $i
  g++ -I /home/wwang/shared/RCRdaemon_Wei_09062013/ -lm $i /home/wwang/shared/RCRdaemon_Wei_09062013/energyStatDaemon.o /home/wwang/shared/RCRdaemon_Wei_09062013/RCR.bb.o /home/wwang/shared/RCRdaemon_Wei_09062013/wrmsr.o  -fopenmp -O3 -o Selected-$i.exe
done
