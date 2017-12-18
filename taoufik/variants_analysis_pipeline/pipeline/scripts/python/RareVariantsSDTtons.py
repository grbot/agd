##Load libs
import getopt
import sys
import argparse
import os
import gzip
import time
import shutil
import re
import fastcsv
##Parser
parser = argparse.ArgumentParser(description='Statistics on Doubletons, Singletons and Tripletons Variants')
parser.add_argument('--freq-input', type=str, default=None, help='Freq file input')
parser.add_argument('--output-folder', type=str, default=None, help='Output folder ')
parser.add_argument('--out', type=str, default=None, help='Output CSV File')
##Parsing inputs
args = parser.parse_args()
outputfolder = args.output_folder
freqinput=args.freq_input
csvoutch = args.out
##Stats
c = 0 #counter of total number of polymorphisms
csingle = 0 #counter of singletons
cdouble = 0 #counter of doubletons
ctriple = 0 #counter of tripletons
spos = []
dpos = []
tpos = []
chrom = re.findall('\\d+', freqinput)[0]
##Folder output
base=os.path.basename(freqinput)



if not outputfolder.endswith("/"):
    outputfolder=outputfolder+"/"
    
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

if not os.path.exists('%s%s/'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0])):
    os.makedirs('%s%s/'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0]))
outFName ='%s%s/%s_%s_Singletons_Doubletons_Tripletons.csv'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],re.findall('\\d+', freqinput)[0])
outFNameReduced ='%s%s/%s_%s_Singletons_Doubletons_Tripletons_Reduced.csv'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],re.findall('\\d+', freqinput)[0])
population = os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0]
##Filtering
with fastcsv.Reader(open(freqinput)) as reader:
    next(reader)
    for row in reader:
        c = c + 1
        if int(row[0].split('\t')[5].split(":")[1]) == 1:
            csingle = csingle + 1
            #print(row)
            spos.append(row[0].split('\t')[1])
        if int(row[0].split('\t')[5].split(":")[1]) == 2:
            cdouble = cdouble + 1
            dpos.append(row[0].split('\t')[1])
        if int(row[0].split('\t')[5].split(":")[1]) == 3:
            ctriple = ctriple + 1
            tpos.append(row[0].split('\t')[1])
 
props = float(csingle)/float(c)
propd = float(cdouble)/float(c)
propt = float(ctriple)/float(c)
outFile = open(csvoutch, 'w')
outFile.write('CHR\tPOS\tSingleton/Doubleton/Tripleton\tTotalSingletons\tPropSingletons\tTotalDoubletons\tPropDoubletons\tTotalTripletons\tPropTripletons\tTotalVariants\n')
for i in range(len(spos)):
	outFile.write('%s\t%s\t%s\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n' % (chrom,spos[i],"S",csingle,props,cdouble,propd,ctriple,propt,c))
for i in range(len(dpos)):
	outFile.write('%s\t%s\t%s\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n' % (chrom,dpos[i],"D",csingle,propd,cdouble,propd,ctriple,propt,c))    
for i in range(len(tpos)):
	outFile.write('%s\t%s\t%s\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n' % (chrom,tpos[i],"T",csingle,props,cdouble,propd,ctriple,propt,c))
outFile.close()
outFile = open(outFNameReduced, 'w')
outFile.write('POP\tCHR\tTotalSingletons\tPropSingletons\tTotalDoubletons\tPropDoubletons\tTotalTripletons\tPropTripletons\tTotalVariants\n')
outFile.write('%s\t%s\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n' % (population,chrom,csingle,props,cdouble,propd,ctriple,propt,c))
outFile.close()
shutil.copy(csvoutch,outFName)

