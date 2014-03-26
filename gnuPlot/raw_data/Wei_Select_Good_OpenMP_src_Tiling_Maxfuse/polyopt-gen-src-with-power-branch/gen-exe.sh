#!/bin/bash
 g++ -I /home/wwang/shared/RCRdaemon_Wei_09062013/ -lm $1 /home/wwang/shared/RCRdaemon_Wei_09062013/energyStatDaemon.o /home/wwang/shared/RCRdaemon_Wei_09062013/RCR.bb.o /home/wwang/shared/RCRdaemon_Wei_09062013/wrmsr.o  -fopenmp -O3
