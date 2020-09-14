###read in the ancestral state of snps
###read in included individuals
###write snp file for importing into adegenet
###edited 12/23/2014, fixed bug on population correlation, loci counting###
###edited 1/9/2015, add a filter on read depth, i.e., only more than certain depth of reads for the loci, the genotype is considered valid ###
###command example: python sampleDownGeno2SFS_NoChrom_Anna.py <popmap> 4,5 <output_filename> <vcf-file> 5
import gzip,sys,random
selspFile = file(sys.argv[1]) #sample to population map
SampleNo = [float(x) for x in sys.argv[2].split(",")] ##down grade to which sample size for the SFS spectrum
prefix = sys.argv[3] #output prefix
infileName = sys.argv[4] #vcf file name  
minDP = int(sys.argv[5]) #minimum read depth for a genotype to be considered
#ancestFile.readline()
#chromAncesAllele = {}
snpLociDict = {}
lociList = []
AlleleStats = {}
#ChromList = sys.argv[7].split(",") #type included chromosome num, separated by ","
#Range = [int(x) for x in sys.argv[8].split(",")] #type range on the chromosome to include, separated by ","
'''
for x in ChromList:
	chromAncesAllele[x]={}
	snpLociDict[x]={}
	lociList[x]=[]
	AlleleStats[x]={}
	

for line in ancestFile:
	(chrom,pos,ancSnp) = line.split()
	pos = int(pos)
	if chrom in chromAncesAllele and pos>=Range[0] and pos <= Range[1] and ancSnp != 'NA':
		chromAncesAllele[chrom][pos] = ancSnp
ancestFile.close()
'''

spPop = {}
for line in selspFile:
	(sp,popName) = line.split()
	spPop[sp] = popName
selspFile.close()
popSet = list(set(spPop.values()))
popSet.sort()
totalSpNo = float(len(spPop))
print "total number of sp is", totalSpNo

outfile = open(prefix + "_MSFS.obs","w")
#outfile.write("//command is: python %s\n"%' '.join(sys.argv))
outfile.write("1 observations. No. of demes and sample sizes are on next line\n")

selsp = [] #record positions of selected individuals


print "convert vcffile " + infileName + " to SFS file for fastsimcoal\n"
if '.gz' in infileName:
	infile = gzip.open(infileName,'r')
else:
	infile = open(infileName,'r')

for line in infile:
	if not line.startswith("##"):
		if line.startswith("#"):
			if not selsp:
				infoList = line.split()
				IndName = infoList[9:]
				for i in xrange(len(IndName)):
					if IndName[i] in spPop:
						#indAllele[IndName[i]] = ''
						selsp.append(i)
						#population.append(spPop[IndName[i]])
			break

#print selsp

for line in infile:
	snpIndList = line.split()
	chrom = snpIndList[0]
	pos = int(snpIndList[1])
	loci = snpIndList[2]
	try:
		snpLociDict[loci].append(pos)
	except:
		snpLociDict[loci] = [pos]
	AlleleStats[pos] = {}
	for x in popSet:
		AlleleStats[pos][x] = [0,0,0,0] #count of 0/0, 0/1, 1/1, and total	
	for x in xrange(9,len(infoList)):
		if x-9 in selsp:
			genotype,dp,ad,lh = snpIndList[x].split(":")
			#genotype = snpIndList[x].split(":")[0]
			if "." not in genotype and int(dp) >= minDP:
				AlleleStats[pos][spPop[IndName[x-9]]][genotype.count("1")] += 1
	for y in popSet:
		AlleleStats[pos][y][3] = sum(AlleleStats[pos][y][:3])
infile.close()
print "done snp calculation"
#print SampleNo
print "popset is ",popSet
multiDimSFS = {}
varLoc = 0
temp = open(prefix + "_pos_count.txt","w")
temp.write("Locus\t" + "\t".join(popSet) + "\n")
#print len(snpLociDict)
for locus in snpLociDict:
	i = 0
	while i< len(snpLociDict[locus]):
		#print "locus is",locus
		#print "snpLociDict[locus] is ",snpLociDict[locus]
		#selNo = random.randint(0,len(snpLociDict[chrom][locus])-1)
		pos = snpLociDict[locus][i]
		allStateValues = [AlleleStats[pos][popSet[j]] for j in xrange(len(popSet))]
		#print allStateValues
		if all(allStateValues[x][3]>=SampleNo[x] for x in xrange(len(allStateValues))): #test whether all pops have adequate samples for the loci
			tempCount = []
			for y in xrange(len(popSet)):
				cat0,cat1,cat2,tt = AlleleStats[pos][popSet[y]]
				tempCount.append(int(round((cat1+cat2*2)*SampleNo[y]/tt)))
			if all(x==0 for x in tempCount) or all(tempCount[x]==int(SampleNo[x]*2) for x in xrange(len(tempCount))):
				i+=1
				#print i
			else:
				#print varLoc
				varLoc += 1
				if sum(tempCount)>sum(SampleNo):
					for x in xrange(len(tempCount)):
						ttt = tempCount[x]
						tempCount[x] = SampleNo[x]*2 - ttt
				temp.write("%s\t%s\n"%(locus,"\t".join([str(x) for x in tempCount])))
				tempCount=tuple(tempCount)
				if tempCount in multiDimSFS:
					multiDimSFS[tempCount]+= 1
				else:
					multiDimSFS[tempCount] = 1
				break		
		else:
			break

outfile.write("%d\t%s\n"%(len(popSet),"\t".join([str(int(x*2)) for x in SampleNo])))
###generate multi-dimensional SFS###
cutTuple = [0]*len(SampleNo)
cutTuple[-1] = -1
maxIt = 1
for x in SampleNo:
	maxIt *= int(x*2+1)

print "total number of MultiDemensional SFS column is ", maxIt
y = 0
while y<maxIt:
	i = -1
	while abs(i)<=len(SampleNo):
		if cutTuple[i]<SampleNo[i]*2:
			cutTuple[i]+=1
			if i<-1:
				for j in xrange(i+1,0):
					cutTuple[j]=0
			break
		i-=1
	#print cutTuple
	x = tuple(cutTuple)
	if x in multiDimSFS:
		outfile.write("%d\t"%multiDimSFS[x])
	else:
		outfile.write("0\t")
	y +=1

outfile.close()
temp.close()

print "total number of variable SNP is %d"%(varLoc)
print "number of invariable loci is %d"%(len(snpLociDict)- varLoc)


