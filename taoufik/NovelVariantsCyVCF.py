##Load libs
import getopt
import sys
##import vcf
from cyvcf2 import VCF, Writer
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
dictionarygnomAD= {'gnomAD_AF','gnomAD_AFR_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_AMR_AF','gnomAD_EAS_AF'}
dictionaryknownsnps={'dbSNPBuildID','CDA','CLNACC','CLNSRC','KGPhase1','KGPilot123','KGPROD','KGValidated','OM','PM','PMC','RS'}
##Print
print 'Filtering %s by %s and gnomAD threshold %f' % (vcfinput,dictionaryknownsnps,thresholdgnomad)
print 'gnomAD dictionary ',dictionarygnomAD,'\n'
print 'filter dictionary ',dictionaryknownsnps,'\n'
##Stats
c = 0 #counter of total number of polymorphisms
novel = 0 #counter of novel variants
afSNPrsID = []
afSNPnovel = []
novelway="Not known snp or threshold gnomad"
##Check novel
def isNovel(record):
	global novelway 
	if record.INFO.get('AF') <= thresholdaf:
		return False
	for val in dictionaryknownsnps:
		if record.INFO.get(val) != None:
			return False
	for val in dictionarygnomAD:
		if record.INFO.get(val) > thresholdgnomad:
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
vcfReader = VCF(vcfinput)
vcfReader.add_info_to_header({'ID': 'NovelSNP', 'Description': 'Detection Way',
    'Type':'Character', 'Number': '1'})
base=os.path.basename(vcfinput)
if os.path.splitext(base)[1]=='.gz':
	outputvcf ='%s%s/%s_Novel.vcf'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
	outputvcfgz ='%s%s/%s_Novel.vcf.gz'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
	outFName = '%s%s/%s_Novel.csv'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
else:
	outputvcf ='%s%s/%s_Novel%s'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(base)[0],os.path.splitext(base)[1])
	outFName = '%s/%s%s_Novel.csv'%(outputfolder,vcfinput.split('_')[0],os.path.splitext(base)[0])
vcfWriter = Writer(outputvcf, vcfReader)
outFile = open(outFName, 'w')
outFile.write('CHR\tPOS\tFREQ\n')
##Filtering
for variant in vcfReader:
	c = c + 1
	#if c==100:
	#	break
	if isNovel(variant): # Variant not in dbSNP and therefore, novel
		variant.INFO["NovelSNP"] = novelway
		vcfWriter.write_record(variant)
		novel = novel + 1
		outFile.write('%s\t%d\t%d\t%f\n' % (variant.CHROM, variant.POS-1, variant.POS, variant.INFO['AF']) )
		afSNPnovel.append(float(variant.INFO['AF']))
	else:
		afSNPrsID.append(float(variant.INFO['AF']))
if os.path.splitext(base)[1]=='.gz':
	with open(outputvcf) as src, gzip.open(outputvcfgz, 'wb') as dst:
		dst.writelines(src)
outFile.close()
print("--- %s seconds ---" % (time.time() - start_time))
sys.stderr.write('Found %d novel variants out of %d (percent: %f)\n' % (novel, c, float(novel)/c ))
sys.stderr.write('Average frequency novel variants: %f\n' % (sum(afSNPnovel)/len(afSNPnovel)) )
sys.stderr.write('Average frequency of known variants: %f\n' % (sum(afSNPrsID)/len(afSNPrsID)) ) 
