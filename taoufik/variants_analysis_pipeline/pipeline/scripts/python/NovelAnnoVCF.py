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
#parser.add_argument('--db-snpid', type=str, default='dbSNPBuildID', help='Filter Not in dbsnp ')
parser.add_argument('--gnom-ad', type=float, default=0.0, help='Filter Not in gnomAD or elsewhere')
parser.add_argument('--af-min', type=float, default=0.0, help='Filter AF ')
parser.add_argument('--output-folder', type=str, default=None, help='Output folder ')
parser.add_argument('--out', type=str, default=None, help='Output vcf file')
##Parsing inputs
args = parser.parse_args()
thresholdaf = args.af_min
thresholdgnomad = args.gnom_ad
outputfolder = args.output_folder
vcfinput=args.vcf_input
vcfoutch = args.out
dictionarygnomAD= {'gnomAD_AF','gnomAD_AFR_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','AF_EXAC','ExAC_AF','ExAC_AFR_AF','ExAC_FIN_AF','ExAC_NFE_AF','ExAC_AMR_AF','ExAC_EAS_AF','AGVP_AF','SAHGP_AF','TRYPANOGEN_AF','KG_AF','KG_AFR_AF','KG_EUR_AF','KG_AMR_AF','KG_EAS_AF'}
dictionaryknownsnps={'dbSNPBuildID','CDA','CLNACC','CLNVI','CLNSIGINCL','CLNSRC','KGPhase1','KGPilot123','KGPROD','KGValidated','OM','PM','PMC','RS','GWASCAT_TRAIT','GWASCAT_REPORTED_GENE','GWASCAT_PUBMED_ID'}

##Print
#print 'Filtering %s by %s and gnomAD threshold %f' % (vcfinput,dictionaryknownsnps,thresholdgnomad)
#print 'gnomAD dictionary ',dictionarygnomAD,'\n'
#print 'filter dictionary ',dictionaryknownsnps,'\n'
##Stats
novelway="Novel"
##Check novel
def isNovel(record):
	global novelway 
	if record.INFO.get('AF') < thresholdaf:
		novelway="Not Novel for having AF < %f" % thresholdaf
		return False
	for val in dictionaryknownsnps:
		if record.INFO.get(val) != None:
			novelway="Not Novel for having  %s" % val
			return False
	for val in dictionarygnomAD:
		if record.INFO.get(val) > thresholdgnomad:
			novelway="Not Novel for having %s > %f" % (val,thresholdgnomad)
			return False
	#if not all(record.INFO[dictionarygnomAD[value]] > self.thresholdgnomad  for value in dictionarygnomAD.keys()):
	return True
##Folder output
base=os.path.basename(vcfinput)
if not outputfolder.endswith("/"):
	outputfolder=outputfolder+"/"	
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
if not os.path.exists('%s%s/'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0])):
	    os.makedirs('%s%s/'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0]))

##Reader and writer
start_time = time.time()
vcfReader = VCF(vcfinput)
vcfReader.add_info_to_header({'ID': 'NovelSNP', 'Description': 'Detection Way',
    'Type':'Character', 'Number': '1'})
if os.path.splitext(base)[1]=='.gz':
	outputvcf ='%s%s/%s_Novel.vcf'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])
	outputvcfgz ='%s%s/%s_Novel.vcf.gz'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0])

else:
	outputvcf ='%s%s/%s_Novel%s'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(base)[0],os.path.splitext(base)[1])

vcfWriter = Writer(outputvcf, vcfReader)
##Filtering
for variant in vcfReader:

	if isNovel(variant):
		variant.INFO["NovelSNP"] = "Novel"
	else:
		variant.INFO["NovelSNP"] = novelway
	vcfWriter.write_record(variant)
	
vcfWriter.close()
if os.path.splitext(base)[1]=='.gz':
	with open(outputvcf) as src, gzip.open(vcfoutch, 'wb') as dst:
		dst.writelines(src)
	os.remove(outputvcf)
	shutil.copy(vcfoutch,outputvcfgz)
#print("--- %s seconds ---" % (time.time() - start_time))
