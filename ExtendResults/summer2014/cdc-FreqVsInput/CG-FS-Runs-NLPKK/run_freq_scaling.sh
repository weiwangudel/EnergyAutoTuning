#!/bin/bash
exe_dir=/home/wwang/shared/grappolo-sb-version-Allan/CG-FS-Runs-NLPKK/

#if [ $# -ne 1 ] ; then
#  echo "provide name!"
#  exit 128 
#fi

module load gcc/4.8.1
./resetDutyCycle.exe
./resetPowerLimit.exe

#for freq in 2701000 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000
for freq in 2700000 2600000 2500000 2400000 2300000 2200000 2100000 1900000 1800000 1700000 1600000 1500000 1400000 1300000 1200000
#for freq in 2700000 2400000 2100000 1800000 1500000 1200000
do 
  echo userspace | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor   
  echo 2700000 | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed

 ##heating CPU
 #export OMP_NUM_THREADS=16
 #for i in `seq 1 1 2` ; do ./covariance_--pluto-ufactor_8_--pluto-fuse_maxfuse_--pluto-parallel_--pluto-tile_32x16x1.exe ; done

 echo userspace | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor   
 echo $freq | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed

 for k in `seq 1 1 1` ; do ./covariance_--pluto-ufactor_8_--pluto-fuse_maxfuse_--pluto-parallel_--pluto-tile_32x16x1.exe ; done

 for num in 16
 do  
   export OMP_NUM_THREADS=$num

     #run five times
     for j in `seq 1 1 5` ; do
     ./graphClustering -f ../grappolo-sb-version-Allan-Wei-Modify-for-OpenMPRuntime/wei_nlpkk.bin -b -c  >> $exe_dir/freq-$freq-num$num.output
     done
   done 

  echo userspace | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor   
  echo 2700000 | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed
done #freq

