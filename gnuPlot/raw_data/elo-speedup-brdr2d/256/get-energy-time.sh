#/bin/bash

grep line * |awk '(NR <=4) {sum1+=$6; sum2+=$10;} (NR>4) {sum3+=$6; sum4+=$10;} (NR==4) {print sum1,sum2} (NR==8) {print sum3,sum4 "\n" sum3/sum1, sum4/sum2}' 
