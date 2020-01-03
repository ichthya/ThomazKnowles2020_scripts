##calculate 95% confidence interval
import sys,math
import numpy as np
prefix = sys.argv[1]
valList = np.genfromtxt(prefix + "_bestlhoods.summary",names = True,delimiter = "\t")
fieldNames = valList.dtype.names
count = int(valList.shape[0])
startPos = int(math.ceil(count*0.025)-1)
endPos = int(math.ceil(count*0.975)-1)
outfile = open(prefix + "_95perc_CI.txt","w")
outfile.write("Parameters\tLowerBound025\tUpperBound975\n")
for y in fieldNames:
	x = valList[y]
	x.sort()
	outfile.write("%s\t%f\t%f\n"%(y,x[startPos],x[endPos]))
outfile.close()
	