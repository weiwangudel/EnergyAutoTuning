#!/bin/bash
exe_dir=`pwd`
#/hpc/shared/home/wwang/mantevo/HPCCG-1.0
#! fix heating issue! Re-run first 50

#if [ $# -ne 1 ] ; then
#  echo "provide name!"
#  exit 128 
#fi


./resetPowerLimit.exe
./resetDutyCycle.exe
./reset*.sh

module load intelc/14.0.2


for i in lu
do 


for dc in `seq 16 -1 6`
do 
 echo "Duty Cycle Value: $dc"
 echo "Compile SetDutyCycle.exe" 
 sed -e  "s/VALUE/"$dc"/g" setDutyCycle-main-template.cpp > setDutyCycle-main.cpp
 g++ -O3 -I/home/wwang/shared/RCRdaemon_Wei_08022014 setDutyCycle-main.cpp  /home/wwang/shared/RCRdaemon_Wei_08022014/libdutyCycle.a -lrt
 mv a.out setDutyCycle.exe 


#for file in `cat covariance-names-tileonly |head -n 9 |tail -n 1`
#do
#  file_name=${file%.*}
#  i=1
#  while [ $i -le 3 ]
#  do
#    echo "Executing $file_name.exe for $i/3 times!"
#    $exe_dir/$file_name.exe 
#    (( i++ ))
#  done
#
#done

./setDutyCycle.exe
 echo "executing during DutyCycle $dc"

     for j in `seq 1 1 10` ; do
        timeout  80	./lulesh.exe 160
        timeout 1800   ./$i-energy.exe   >> $exe_dir/$i-dc-$dc.output
#      timeout 4 ./lulesh.exe 120     # clear cache to avoid as much as variance
#      ./test_HPCCG 200 200 200 >> $exe_dir/dc-200-$dc.output	
     done

./resetDutyCycle.exe

done 

done   #each benchmark
