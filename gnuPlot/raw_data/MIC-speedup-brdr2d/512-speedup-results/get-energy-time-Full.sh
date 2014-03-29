#/bin/bash

grep App *.output* |awk '(NR <=6) {sum1+=$4; sum2+=$8;} (NR>6) {sum3+=$4; sum4+=$8;} (NR==6) {print sum1,sum2} (NR==12) {print sum3,sum4 "\n" sum3/sum1, sum4/sum2}' 
