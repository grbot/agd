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
	row = ['CHR\tPOS\tSingleton/Doubleton/Tripleton\tTotalSingletons\tPropSingletons\tTotalDoubletons\tPropDoubletons\tTotalTripletons\tPropTripletons\tTotalVariants']	
		
	with fastcsv.Writer(open(currentpop+pop+"_All.csv", 'w')) as writer:
		writer.writerow(row)
		for freq in popfreqs:
			with fastcsv.Reader(open(currentpop+freq)) as reader:
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
	row = ['POP\tCHR\tTotalSingletons\tPropSingletons\tTotalDoubletons\tPropDoubletons\tTotalTripletons\tPropTripletons\tTotalVariants']	
		
	with fastcsv.Writer(open(currentpop+pop+"_All_Reduced.csv", 'w')) as writer:
		writer.writerow(row)
		for freq in popfreqs:
			with fastcsv.Reader(open(currentpop+freq)) as reader:
				next(reader)
				for row in reader:
					writer.writerow(row)

	


