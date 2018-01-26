##Load libs
import getopt
import sys
from cyvcf2 import VCF, Writer
import argparse
import os
import gzip
import time
import shutil
##Parser
parser = argparse.ArgumentParser(description='Filtering Novel SNPs')
parser.add_argument('--vcf-input', type=str, default=None, help='Output Filtred')
parser.add_argument('--output-folder', type=str, default=None, help='Output folder ')
parser.add_argument('--out', type=str, default=None, help='Output CSV file')
##Parsing inputs
args = parser.parse_args()
outputfolder = args.output_folder
vcfinput=args.vcf_input
csvoutch = args.out
##Stats
c = 0 #counter of total number of polymorphisms
novel = 0 #counter of novel variants
afSNPrsID = []
afSNPnovel = []
SNPnovelchrom = []
SNPnovelpos = []
afSNPnovelpass5percent = []
afSNPnovelpass10percent = []
afSNPnovelprop5percent = []
afSNPnovelprop10percent = []
c5=0
c10=0
##Folder output
base=os.path.basename(vcfinput)
if not outputfolder.endswith("/"):
	outputfolder=outputfolder+"/"	
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
if not os.path.exists('%s%s/'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0])):
	    os.makedirs('%s%s/'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0]))

##Reader and writer
vcfReader = VCF(vcfinput)
if os.path.splitext(base)[1]=='.gz':
	outFName ='%s%s/%s_Stats.csv'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
else:
	outFName ='%s%s/%s_Stats.csv'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(base)[0],os.path.splitext(base)[1])
##Filtering
for variant in vcfReader:
    c = c+1
    #print variant.INFO["NovelSNP"]
    if variant.INFO["NovelSNP"]=="Novel":
    	novel = novel +1
    	afSNPnovel.append(variant.INFO['AF'])
    	SNPnovelchrom.append(variant.CHROM)
    	SNPnovelpos.append(variant.POS)
    	if variant.INFO['AF'] > 0.05:
    		afSNPnovelpass5percent.append("PASS")
    		c5 = c5+1
    	else:
    		afSNPnovelpass5percent.append("REJECT")
    	if variant.INFO['AF'] > 0.1:
    		afSNPnovelpass10percent.append("PASS")
    		c10 = c10+1
    	else:
    		afSNPnovelpass10percent.append("REJECT")
    else:
    	print variant.INFO["NovelSNP"]
    	afSNPrsID.append(variant.INFO['AF'])
outFile = open(csvoutch, 'w')
outFile.write('CHR\tPOS\tFREQ\tPROP5\tPASS5\tPROP10\tPASS10\tAVFREQNOVEL\tAVFREQKNOWN\tRATIONOVEL\tNUMNOVEL\tNUMKNOWN\tNUMTOTAL\n')
if novel == 0:
	propnovel = 0
	avgfreqnovel = 0
	avgfreqknown = sum(afSNPrsID)/len(afSNPrsID)
	prop5 = 0
	prop10 = 0
	cknown = c - novel
	outFile.write('%s\t%d\t%f\t%f\t%s\t%f\t%s\t%f\t%f\t%f\t%f\t%d\t%d\n' %("-",0,0,prop5,"-",prop10,"-",avgfreqnovel,avgfreqknown,propnovel, novel,cknown,c))
else:
	if c==novel:
		propnovel = float(novel)/float(c)
		avgfreqnovel = sum(afSNPnovel)/len(afSNPnovel)
		avgfreqknown = 0
		prop5 = float(c5)/float(novel)
		prop10 = float(c10)/float(novel)
		cknown = c - novel
	else:
		propnovel = float(novel)/float(c)
		avgfreqnovel = sum(afSNPnovel)/len(afSNPnovel)
		avgfreqknown = sum(afSNPrsID)/len(afSNPrsID)
		prop5 = float(c5)/float(novel)
		prop10 = float(c10)/float(novel)
		cknown = c - novel
	for i in range(len(afSNPnovel)):
		outFile.write('%s\t%d\t%f\t%f\t%s\t%f\t%s\t%f\t%f\t%f\t%f\t%d\t%d\n' % (SNPnovelchrom[i],SNPnovelpos[i],afSNPnovel[i],prop5,afSNPnovelpass5percent[i],prop10,afSNPnovelpass10percent[i],avgfreqnovel,avgfreqknown,propnovel, novel,cknown,c))
##Writing Stats
outFile.close()

shutil.copy(csvoutch,outFName)

