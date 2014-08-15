


for num in 16 15 8 
do 
   for file in fuse_nofuse_parallel_tile_32x16x1.data
   do
    
    awk '($3=='$num') {print $2 " " $3 " " $5 " " $5*$4 " " $5*$4*$4}'  $file > temp
    
    minE=$(sort -k 3 -r -n temp |tail -n 1 | cut -d" " -f3)
    minED=$(sort -k 4 -r -n temp |tail -n 1 | cut -d" " -f4)
    minEDD=$(sort -k 5 -r -n temp |tail -n 1 | cut -d" " -f5)
    
   awk '{print $1 " " $2 " " $3/'$minE' " " $4/'$minED' " " $5/'$minEDD' }' temp > gnuplot.data 
   tac gnuplot.data > temp
   mv temp gnuplot.data
   gnuplot draw.gnpl
   mv output.pdf $num.pdf
   done 
done 
