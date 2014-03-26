#!/bin/bash

for i in `ls max*.c `
do
#  sed -i 's/Real_t tx;/energyDaemonInit();\nReal_t tx;/g' $i
#  sed -i 's/return 0;/energyDaemonTerm();\n return 0;/g' $i
#  sed -i 's/Index_t jj;/Index_t jj;\n energyDaemonEnter();/g' $i
#  sed -i 's/if (hgcoef/energyDaemonExit("file", 1);\n if (hgcoef/g' $i
#  sed -i 's/1.0e-6/1.0e-4/g' $i
#  g++ -I /home/wwang/shared/RCRdaemon_Wei_09062013/ -lm $i /home/wwang/shared/RCRdaemon_Wei_09062013/energyStatDaemon.o /home/wwang/shared/RCRdaemon_Wei_09062013/RCR.bb.o /home/wwang/shared/RCRdaemon_Wei_09062013/wrmsr.o  -fopenmp -O3 -o Selected-$i.exe
  g++ -I /home/wwang/shared/RCRdaemon_Wei_09062013/ -lm $i /home/wwang/shared/RCRdaemon_Wei_09062013/energyStatDaemon.o /home/wwang/shared/RCRdaemon_Wei_09062013/RCR.bb.o -fopenmp -O3 -o Selected-$i.exe
done
