


for num in `seq 16 -1 16`
do 
   for file in data
   do
    
    awk '($2=='$num') {print $1 " " $2 " " $4 " " $4*$3 " " }'  $file > temp
    
    minE=$(sort -k 3 -r -n temp |tail -n 1 | cut -d" " -f3)
    minED=$(sort -k 4 -r -n temp |tail -n 1 | cut -d" " -f4)
    
    
   awk '{print $1 " " $2 " " $3/'$minE' " " $4/'$minED' }' temp > gnuplot.data 
   tac gnuplot.data > temp
   mv temp gnuplot.data
   gnuplot draw.gnpl
   mv output.pdf $num.pdf
   done 
done 
