##Load libs
import getopt
import sys
from cyvcf2 import VCF, Writer
import argparse
import os
import gzip
import time
import shutil
import csv
import re
##Parser
parser = argparse.ArgumentParser(description='Rare Variants Impact')
parser.add_argument('--vcf-input', type=str, default=None, help='VCF input')
parser.add_argument('--rarevariants-input', type=str, default=None, help='Freq Count input ')
parser.add_argument('--output-folder', type=str, default=None, help='Output Folder ')
##Parsing inputs
args = parser.parse_args()
outputfolder = args.output_folder
vcfinput=args.vcf_input
rarevariantsinput=args.rarevariants_input
ontoeffects={'coding_sequence_variant','chromosome','duplication','inversion','coding_sequence_variant','inframe_insertion','disruptive_inframe_insertion','inframe_deletion','disruptive_inframe_deletion','downstream_gene_variant','exon_variant','exon_loss_variant','exon_loss_variant','duplication','duplication','inversion','inversion','frameshift_variant','gene_variant','feature_ablation','duplication','gene_fusion','gene_fusion','bidirectional_gene_fusion','rearranged_at_DNA_level','intergenic_region','conserved_intergenic_variant','intragenic_variant','intron_variant','conserved_intron_variant','miRNA','missense_variant','initiator_codon_variant','stop_retained_variant','protein_protein_contact','structural_interaction_variant','rare_amino_acid_variant','splice_acceptor_variant','splice_donor_variant','splice_region_variant','splice_region_variant','splice_region_variant','stop_lost','5_prime_UTR_premature','start_codon_gain_variant','start_lost','stop_gained','synonymous_variant','start_retained','stop_retained_variant','transcript_variant','feature_ablation','regulatory_region_variant','upstream_gene_variant','3_prime_UTR_variant','3_prime_UTR_truncation','exon_loss','5_prime_UTR_variant','5_prime_UTR_truncation','exon_loss_variant','sequence_feature','exon_loss_variant'}
chrom = re.findall('\\d+', vcfinput)[0]
##Folder output
base=os.path.basename(vcfinput)
if not outputfolder.endswith("/"):
    outputfolder=outputfolder+"/"
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
outFName ='%s%s/%s_%s_Effect_Impact.csv'%(outputfolder,os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],os.path.splitext(os.path.splitext(base)[0])[0].split('_')[0],chrom)
print outFName
def formatANN(record):
	formattedann =""
	if record.INFO.get("LOF")==None:
		listeffects= record.INFO.get("ANN").split(',') 
		for val in ontoeffects:
			effecto = "None"
			for effect in listeffects:
				if effect.split('|')[1]==val:
					effecto = effect.split('|')[2]
					break
			formattedann = formattedann + '\t'+effecto	
	else:
		for val in ontoeffects:
			formattedann = formattedann + '\tLOF'	
			
	#if not all(record.INFO[dictionarygnomAD[value]] > self.thresholdgnomad  for value in dictionarygnomAD.keys()):
	return formattedann
##Reader and writer
vcf = VCF(vcfinput)
header = ""
for effect in ontoeffects:
	header = header + "\t"+effect
header = "CHROM\tPOS"+header+"\n"
##Filtering
c = 0
with open(rarevariantsinput) as rcsvfile:
	reader = csv.reader(rcsvfile, delimiter='\t')
	next(reader)
	outfn = open(outFName,'w')
	outfn.write(header)
	for row in reader:
		vio = vcf("%d:%d-%d"%(int(row[0]),int(row[1]),int(row[1])))
		for v in vio:
			outfn.write("%d\t%d%s\n"%(int(row[0]),int(row[1]),formatANN(v)))
	outfn.close()
