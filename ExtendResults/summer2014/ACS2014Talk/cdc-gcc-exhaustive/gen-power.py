import sys
# Given a list that contains N items, print 2**N items out of them
def printPowerofList (listA):
#    print listA
#    for a in listA:
#        print a,
    for i in range(pow(2, len(listA))):
       intList = [int(x) for x in bin(i)[2:].zfill(len(listA))]
       #print intList
       flagList = []
       for j in range (len(intList)):
           if (intList[j] == 1):
                flagList += listA[j].rstrip(), ## get rid of \n
       for k in flagList:
           print k,
       print '\n',
         
filename = sys.argv[1]       
listA = []
# Echo the contents of a file
f = open(filename, 'rU')
for line in f:   ## iterates over the lines of the file
  listA += line,    ## trailing , so print does not add an end-of-line char
                 ## since 'line' already includes the end-of line.
f.close()
printPowerofList(listA)

    
