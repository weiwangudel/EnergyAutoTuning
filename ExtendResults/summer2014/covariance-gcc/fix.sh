#!/bin/bash

usage () {
	echo "Usage: $0 optimization name (e.g. maxfuse fuse 2mm)"
        exit 1
}

[[ $# -eq 0 ]] && usage



grep $1 data  > with-$1-1


sed -i "s/$3_--pluto-$2_$1_--pluto-tile//g" with-$1-1

grep $3 with-$1-1 > with-$1-notile

sed -i "s/$3_--pluto-$2_$1.exe//g" with-$1-notile

grep $3 with-$1-notile > with-$1

############# Over getting with, begin getting without #############

cut -d" " -f1 with-$1 > temp

sed -i "s/_--pluto-$2_$1//g" temp

mv temp wo-$1-names

echo "begin greping performance without $1..."
for i in `cat wo-$1-names`
do 
 grep $i data >> wo-$1
done 

paste with-$1 wo-$1 > with-wo-$1
