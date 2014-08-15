#!/bin/bash
exe_dir=/home/wwang/shared/grappolo-sb-version-Allan/FreqScalingRuns

#if [ $# -ne 1 ] ; then
#  echo "provide name!"
#  exit 128 
#fi

module load gcc/4.8.1

for freq in 2701000 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000
#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
do 
  echo userspace | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor   
  echo $freq | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed

 #heating CPU
 export OMP_NUM_THREADS=16
 for i in `seq 1 1 20` ; do ./covariance_--pluto-ufactor_8_--pluto-fuse_maxfuse_--pluto-parallel_--pluto-tile_32x16x1.exe ; done

 for num in 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
 do  
   export OMP_NUM_THREADS=$num

   for file in graphClustering 
   do
     #run five times
     for j in `seq 1 1 5` ; do
     $exe_dir/$file -f pnnl.bin -b -c  >> $exe_dir/freq-$freq-num$num.output
     done
   done 

  done  #num
done #freq

  echo userspace | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor   
  echo 2700000 | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed
