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
##Parser
parser = argparse.ArgumentParser(description='Combine Statistics on Doubletons, Singletons and Tripletons Variants')
parser.add_argument('--freqs-folder', type=str, default=None, help='Freq files inputs')
parser.add_argument('--output-folder', type=str, default=None, help='Output folder ')
args = parser.parse_args()
freqsfolder = args.freqs_folder
if not freqsfolder.endswith("/"):
	freqsfolder=freqsfolder+"/"
popdirs = os.listdir(freqsfolder)
currentpop = freqsfolder
for pop in popdirs:
	#print "Pop : ",pop
	currentpop = freqsfolder+pop+"/"
	#print "Current Pop folder : ", currentpop
	popfreqs = [f for f in os.listdir(currentpop) if re.match(r'.*Singletons_Doubletons_Tripletons\.csv', f)]
	#popfreqs = os.listdir(currentpop)
	#print popfreqs
	row = ['CHR','POS','Singleton/Doubleton/Tripleton','TotalSingletons','PropSingletons','TotalDoubletons','PropDoubletons','TotalTripletons','PropTripletons','TotalVariants']	
	#row = ['POP','CHR','TotalSingletons','PropSingletons','TotalDoubletons','PropDoubletons','TotalTripletons','PropTripletons','TotalVariants']$		
	with open(currentpop+pop+"_All.csv", 'w') as wcsvfile:
		writer = csv.writer(wcsvfile, delimiter='\t')	
		writer.writerow(row)
		for freq in popfreqs:
			with open(currentpop+freq) as rcsvfile:
				reader = csv.reader(rcsvfile, delimiter='\t')
				next(reader)
				for row in reader:
					writer.writerow(row)
for pop in popdirs:
	#print "Pop : ",pop
	currentpop = freqsfolder+pop+"/"
	#print "Current Pop folder : ", currentpop
	popfreqs = [f for f in os.listdir(currentpop) if re.match(r'.*Singletons_Doubletons_Tripletons_Reduced\.csv', f)]
	#popfreqs = os.listdir(currentpop)
	#print popfreqs
	row = ['POP','CHR','TotalSingletons','PropSingletons','TotalDoubletons','PropDoubletons','TotalTripletons','PropTripletons','TotalVariants']	
	with open(currentpop+pop+"_All_Reduced.csv", 'w') as wcsvfile:
		writer = csv.writer(wcsvfile, delimiter='\t')	
		writer.writerow(row)
		for freq in popfreqs:
			with open(currentpop+freq) as rcsvfile:
				reader = csv.reader(rcsvfile, delimiter='\t')
				next(reader)
				for row in reader:
					writer.writerow(row)

	


