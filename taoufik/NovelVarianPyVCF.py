##Load libs
import getopt
import sys
import vcf
import argparse
import os
import gzip
import time
##Parser
parser = argparse.ArgumentParser(description='Filtering Novel SNPs')
parser.add_argument('--vcf-input', type=str, default=None, help='Output Filtred')
#parser.add_argument('--db-snpid', type=str, default='dbSNPBuildID', help='Filter Not in dbsnp ')
parser.add_argument('--gnom-ad', type=float, default=0.0, help='Filter Not in gnomAD ')
parser.add_argument('--af-min', type=float, default=0.0, help='Filter AF ')
parser.add_argument('--output-folder', type=str, default=None, help='Output folder ')
##Parsing inputs
args = parser.parse_args()
thresholdaf = args.af_min
thresholdgnomad = args.gnom_ad
outputfolder = args.output_folder
vcfinput=args.vcf_input
dictionarygnomAD= {'gnomAD_AF','gnomAD_AFR_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','AF_EXAC','ExAC_AF','ExAC_AFR_AF','ExAC_FIN_AF','ExAC_NFE_AF','ExAC_AMR_AF','ExAC_EAS_AF','AGVP_AF','SAHGP_AF','TRYPANOGEN_AF','KG_AF','KG_AFR_AF','KG_EUR_AF','KG_AMR_AF','KG_EAS_AF'}
dictionaryknownsnps={'dbSNPBuildID','CDA','CLNACC','CLNVI','CLNSIGINCL','CLNSRC','KGPhase1','KGPilot123','KGPROD','KGValidated','OM','PM','PMC','RS','GWASCAT_TRAIT','GWASCAT_REPORTED_GENE','GWASCAT_PUBMED_ID'}
##Print
print 'Filtering %s by %s and gnomAD threshold %f' % (vcfinput,dictionaryknownsnps,thresholdgnomad)
print 'gnomAD dictionary ',dictionarygnomAD,'\n'
print 'filter dictionary ',dictionaryknownsnps,'\n'
##Stats
c = 0 #counter of total number of polymorphisms
novel = 0 #counter of novel variants
afSNPrsID = []
afSNPnovel = []
##Check novel
def isNovel(record):
	recordinfokeys = record.INFO.keys()
	if 'AF' in recordinfokeys:
		if record.INFO['AF'] < thresholdaf:
			return False
	else:
		return False 
	if len(set.intersection(dictionaryknownsnps,record.INFO.keys())) > 1:
		return False
	for val in dictionarygnomAD:
		if val in recordinfokeys:
			if record.INFO[val] > thresholdgnomad:
				return False
	#if not all(record.INFO[dictionarygnomAD[value]] > self.thresholdgnomad  for value in dictionarygnomAD.keys()):
	return True
##Folder output
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
if not os.path.exists('%s%s/'%(outputfolder,vcfinput.split('_')[0])):
    os.makedirs('%s%s/'%(outputfolder,vcfinput.split('_')[0]))

##Reader and writer
start_time = time.time()
vcfReader = vcf.Reader(filename=vcfinput)
base=os.path.basename(vcfinput)
if os.path.splitext(base)[1]=='.gz':
	outputvcf ='%s%s/%s_Novel.vcf'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
	outputvcfgz ='%s%s/%s_Novel.vcf.gz'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
	outFName = '%s%s/%s_Novel.csv'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
else:
	outputvcf ='%s%s/%s_Novel%s'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(base)[0],os.path.splitext(base)[1])
	outFName = '%s/%s%s_Novel.csv'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(base)[0])
vcfWriter = vcf.Writer(file(outputvcf, 'w'), vcfReader)
outFile = open(outFName, 'w')
outFile.write('CHR\tPOS\tFREQ\n')
##Filtering
for variant in vcfReader:
	c = c + 1
	#if c==100:
	#	break
	if isNovel(variant): # Variant not in dbSNP and therefore, novel
		vcfWriter.write_record(variant)
		novel = novel + 1
		outFile.write('%s\t%d\t%d\t%f\n' % (variant.CHROM, variant.POS-1, variant.POS, variant.INFO['AF'][0]) )
		afSNPnovel.append(float(variant.INFO['AF'][0]))
	else:
		afSNPrsID.append(float(variant.INFO['AF'][0]))
if os.path.splitext(base)[1]=='.gz':
	with open(outputvcf) as src, gzip.open(outputvcfgz, 'wb') as dst:
		dst.writelines(src)
outFile.close()
print("--- %s seconds ---" % (time.time() - start_time))
sys.stderr.write('Found %d novel variants out of %d (percent: %f)\n' % (novel, c, float(novel)/c ))
sys.stderr.write('Average frequency novel variants: %f\n' % (sum(afSNPnovel)/len(afSNPnovel)) )
sys.stderr.write('Average frequency of known variants: %f\n' % (sum(afSNPrsID)/len(afSNPrsID)) ) 
