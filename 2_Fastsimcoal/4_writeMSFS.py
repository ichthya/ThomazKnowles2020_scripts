###execute fastsimcoal
import glob,os,sys
import subprocess as sp
prefix = sys.argv[1]
rep = int(sys.argv[2])
obsf = file(prefix + "_sim_MSFS.obs")
obsAll = obsf.readlines()
outfile = open(prefix + "_MSFS.obs","w")
outfile.write("1 observations. No. of demes and sample sizes are on next line\n%s%s"%(obsAll[1],obsAll[rep+1]))
outfile.close()
