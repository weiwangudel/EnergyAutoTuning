import sys
listA = []

# Echo the contents of a file
f = open(sys.argv[1], 'rU')
for line in f:   ## iterates over the lines of the file
  listA += line,    ## trailing , so print does not add an end-of-line char
                 ## since 'line' already includes the end-of line.
f.close()
i = 0;
for each in listA:
#    print bin(i)[2:].zfill(8), 
    print bin(int(each))[2:].zfill(8), 
    print each,                 
    i += 1;

    
