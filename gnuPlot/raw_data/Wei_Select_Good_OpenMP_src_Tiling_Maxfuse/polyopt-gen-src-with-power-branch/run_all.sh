#!/bin/bash

for i in `ls Sele* |head -n 1`
do
  j=1
  while [ $j -le 10 ] 
  do
    ./$i
    (( j++ ))
  done
done

for i in `ls Sele*`
do
./$i > $i.output
done
