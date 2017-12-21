##Load libs
import getopt
import sys
import argparse
import os
import gzip
import time
import shutil
import re
import csv
import subprocess
##Parser
parser = argparse.ArgumentParser(description='Batch Variants Effects Impact')
parser.add_argument('--freqs-folder', type=str, default=None, help='Freq Files Inputs')
parser.add_argument('--vcfs-folder', type=str, default=None, help='Vcf Files Inputs ')
parser.add_argument('--output-folder', type=str, default=None, help='Output Folder ')
parser.add_argument('--python-script', type=str, default=None, help='Variants Effects Python Script ')
args = parser.parse_args()
freqsfolder = args.freqs_folder
vcfsfolder = args.vcfs_folder
#print vcfsfolder
#print freqsfolder
outputfolder = args.output_folder
pythonscript = args.python_script
if not freqsfolder.endswith("/"):
	freqsfolder=freqsfolder+"/"
if not vcfsfolder.endswith("/"):
	vcfsfolder=vcfsfolder+"/"
if not outputfolder.endswith("/"):
	outputfolder=outputfolder+"/"
popdirs = os.listdir(freqsfolder)
currentpopfreq = freqsfolder
currentpopvcf = vcfsfolder
#print vcfsfolder
#print freqsfolder
outsh = open("rarebatch.sh",'w')
outsh.write("#!/bin/bash\n")
for pop in popdirs:
	#print "Pop : ",pop
	currentpopfreq = freqsfolder+pop+"/"
	currentpopvcf = vcfsfolder+pop+"/"
	#print "Current Pop folder : ", currentpop
	popfreqs = [f for f in os.listdir(currentpopfreq) if re.match(r'.*Singletons_Doubletons_Tripletons\.csv', f)]
	popvcfs = [f for f in os.listdir(currentpopvcf) if re.match(r'.*Novel.vcf.gz$', f)]
	popfreqs.sort()
	popvcfs.sort()
	#print currentpopvcf
	#print currentpopfreq
	#print popfreqs
	#print popvcfs
	for i in range(len(popfreqs)):
		#print popvcfs[i],"and",popfreqs[i]
		bashCommand = "echo \"python %s --vcf-input %s --rarevariants-input %s --output-folder %s\" | qsub -N %s -l nodes=1:ppn=1,walltime=2:00:00,mem=5GB -e %s -o %s\n" % (pythonscript,currentpopvcf+popvcfs[i],currentpopfreq+popfreqs[i],outputfolder,"EEFin"+pop,"error-EEFin"+pop+".log","out-EEFin"+pop+".log")
		outsh.write(bashCommand)
		#print bashCommand
		#process = subprocess.Popen(bashCommand.split(),shell = True, stdout=subprocess.PIPE)
		#output, error = process.communicate()
		#print output
outsh.close()



	


