"""
Author: Mamana M.
Affiliation: University of Cape Town
Aim: A simple Snakemake workflow for GAPW african genomic diversity stream.
Date: Mon June 21 14:03:11 CET 2016
Run: snakemake -s afr_genomic_diversity.snakefile --configfile afr_genomic_diversity.yaml --timestamp -p --cores

Latest modification:
  - TODO
"""


import sys, os, subprocess, time, gzip, glob
import fileinput as fi
from snakemake.utils import R

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")

def get_envs():
    return str(subprocess.Popen("conda info --root".split(), stdout=subprocess.PIPE).communicate()[0].strip(), "utf-8")

##-----------------------------------------------##
## Variable from config file                     ##
##                                               ##
##-----------------------------------------------##

basedir = config['basedir']
homedir = config['homedir']
VENV_DIR = get_envs()+"/envs/"+config["venv_dir"]+"/bin"
data_dir = config['data_dir']
genomedata_path = config['genomedata_path']
prefix = config["prefix"]
suffix = config["suffix"]
log_dir = basedir+"/LOG"
sample_file = config["sample_file"]
ref_file_format = config['ref_file_format']
snpEff = config['snpEff']
snpSift = config['snpSift']
snpEff_human_db = config['snpEff_human_db']
snpEff_database = config['snpEff_database']
snpEff_dbsnp_url = config['snpEff_dbsnp_url']
dbsnp_vcf = config['dbsnp_vcf']
snpEff_gwascatalog = config['snpEff_gwascatalog']
POPS = {pop:config['POPS'][pop] for pop in config['POPS']}

POPS_ = {}
for pop in open(sample_file):
    pop = pop.strip().split()[1]
    if pop not in POPS_:
        POPS_[pop] = 1
    else:
        POPS_[pop] += 1
POP_samples = [pop+':'+(str(POPS_[pop])) for pop in POPS_]

message("The current working directory is " + homedir)
message("The data directory is " + data_dir)
try:
    snakemake.utils.makedirs(log_dir)
except:
    pass
finally:
    message("Log stored in " + log_dir)


# CHRS = [str(chr) for chr in list(range(int(config['chromosomes'].split('-')[0]), int(config['chromosomes'].split('-')[-1])+1))]
CHRS = [chrm.strip() for chrm in config['chromosomes'].split(',')]

try:
    chrXY = config['chromosomesXY'].split(',')
except:
    chrXY = []

message("Chromosomes to be used: " + ', '.join(CHRS+chrXY))
message("Populations in sample file: " + ', '.join(POP_samples))
##-----------------------------------------------##
## A set of functions
##-----------------------------------------------##


def ind_pop_id_filter(sample_file_in, sample_file_out):
    """
    """
    lineno = 0
    message("Reading "+sample_file_in)
    message("Writing "+sample_file_out)
    sample_file_out_ = open(sample_file_out, 'w')
    for line in open(sample_file_in):
        if lineno > 1:
            line = line.strip().split()
            sample_file_out_.writelines('\t'.join([line[0], line[2]])+'\n')
        lineno += 1
    sample_file_out_.close()

def generate_dp_report(fil_bef, fil_after, dp_report, CHRS=CHRS):
    """
    """
    print(fil_bef)

def split_sample_list_per_pop(POP_, input_file, output_file):
    '''
    Read integrated_call_samples_v3.20130502.ALL.panel and split it by population
    '''
    i = 1 #line no
    data = []
    message("Extracting "+POP_+' to '+output_file)
    for line in open(input_file):
        if i > 1:
            line = line.strip().split()
            POP = line[1].strip()
            if POP == POP_:
                data.append('\t'.join(line)+'\n')
        i += 1
    if len(data) > 0:
        out = open(output_file, "w")
        for elt in data:
            out.writelines(elt)
        out.close()

rule all:
    '''
    To run everything to the end of the Pipeline
    '''
    input:
        expand(data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.daf.frq", POP=POPS),
        expand(data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.daf.frq", chrm=CHRS, POP=POPS),

rule filter_vcf_by_depth_all:
    '''
    |-- Step 1: Filter by depth
    '''
    input:
        expand(data_dir+"/VCF_FILTERED/baylor_phased_chr{chrm}_dp6.vcf.gz", chrm=CHRS),

rule add_ANC_to_VCF_all:
    '''
    |-- Step 2: Fill ancestral allele using in-house python script
    '''
    input:
        expand(data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc.vcf.gz", chrm=CHRS),
        expand(data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc_f.vcf.gz", chrm=CHRS),
        #with 1000G data using fill-aa of vcftools NOT WORKING :-)
        # expand(data_dir+"/reference/ancestral_alignments/human_ancestor_{chrm}.fa.gz", chrm=CHRS),
        # expand(data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_AA.vcf.gz", chrm=CHRS),

rule annotate_vcf_snpEff_all:
    '''
    |-- Step 3: Annotate VCF using snpEff
    '''
    input:
        snpEff_database+"/"+snpEff_human_db+".txt",
        snpEff_database+"/gwascatalog.txt",
        expand(data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp.vcf.gz", chrm=CHRS),
        expand(data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz", chrm=CHRS),

rule daf_by_pop_all:
    '''
    |-- Step 4: DAF per population and per chrm using VCFTools
    '''
    input:
        expand(data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.daf.frq", POP=POPS),
        expand(data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.daf.frq", chrm=CHRS, POP=POPS),

rule compress_vcf_all:
    '''
    Filter by depth
    '''
    input:
        expand(data_dir+"/VCF/"+prefix+"{chrm}"+suffix+".vcf.gz", chrm=CHRS),
        expand(data_dir+"/VCF/"+prefix+"{chrm}"+suffix+".vcf.gz.tbi", chrm=CHRS),

rule split_POP_samples_all:
    '''
    |-- Fill ancestral allele using in-house python script
    '''
    input:
        expand(basedir+"/samples/{POP}.sample", POP=POPS),

rule split_vcf_by_pop_all:
    '''
    |-- Split vcf per population
    '''
    input:
        expand(data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz", chrm=CHRS, POP=POPS),

rule concat_vcf_in_pop_all:
    '''
    |-- Concat VCF chrm in population using BCFTools
    '''
    input:
        expand(data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.vcf.gz", POP=POPS),

rule ind_pop_id_filter:
    '''
    '''
    input:
        sample_file = sample_file
    output:
        sample_file = os.path.splitext(sample_file)[0]+"_ind_pop_id_only"+os.path.splitext(sample_file)[1]
    params:
        cluster = config['queue_small']+" -o "+basedir+"/LOG/ind_pop_id_filter.out -e "+basedir+"/LOG/ind_pop_id_filter.err -N ind_pop_id_filter"
    run:
        ind_pop_id_filter(input.sample_file, output.sample_file)

rule compress_vcf:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF/"+prefix+"{chrm}"+suffix+".vcf"
    output:
        vcf_file = data_dir+"/VCF/"+prefix+"{chrm}"+suffix+".vcf.gz"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/compress_{chrm}.out -N compress_{chrm}",
    run:
        shell("bgzip -f {input.vcf_file}")

rule index_raw_vcf:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF/"+prefix+"{chrm}"+suffix+".vcf.gz"
    output:
        vcf_file = data_dir+"/VCF/"+prefix+"{chrm}"+suffix+".vcf.gz.tbi"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/tbi_raw_{chrm}.out -N tbi_raw_{chrm}",
    run:
        shell("bcftools index --tbi {input.vcf_file}")

# rule index_vcf:
#     '''
#     '''
#     input:
#         vcf_file = data_dir+'/VCF/'+prefix+"{chrm}"+suffix+".vcf.gz"
#     output:
#         vcf_file = data_dir+'/VCF/'+prefix+"{chrm}"+suffix+".vcf.gz.tbi"
#     params:
#         cluster = config['queue_small']+" -o "+basedir+"/LOG/index_vcf_{chrm}.out -e "+basedir+"/LOG/index_vcf_{chrm}.err -N index_vcf_{chrm}",
#         vcf_file = data_dir+'/VCF/'+prefix+"{chrm}"+suffix+".vcf"
#     run:
#         shell("gunzip -f {input.vcf_file}")
#         shell("bgzip -f {params.vcf_file}")
#         shell("bcftools index --tbi {input.vcf_file}")

rule filter_vcf_by_depth:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF/"+prefix+"{chrm}"+suffix+".vcf.gz"
    output:
        vcf_file = data_dir+"/VCF_FILTERED/baylor_phased_chr{chrm}_dp6.vcf.gz"
    params:
        cluster = config['queue_small']+" -o "+basedir+"/LOG/filter_vcf_by_depth_{chrm}.out -e "+basedir+"/LOG/filter_by_depth_{chrm}.err -N filter_by_depth_{chrm}",
    run:
        shell(VENV_DIR+"/bcftools \
            view -i \'DP>6\' \
            {input.vcf_file}\
            | "+VENV_DIR+"/bgzip -c > {output.vcf_file}")
        # shell(VENV_DIR+"/vcftools \
        #     --gzvcf {input.vcf_file} \
        #     --min-meanDP 6 \
        #     --temp "+data_dir+" \
        #     --recode --stdout | gzip -c > {output.vcf_file}")

rule index_vcf_after_dp:
    '''
    '''
    input:
        vcf_file = data_dir+'/VCF_FILTERED/'+prefix+"{chrm}"+suffix+"_dp6.vcf.gz"
    output:
        vcf_file = data_dir+'/VCF_FILTERED/'+prefix+"{chrm}"+suffix+"_dp6.vcf.gz.tbi"
    params:
        cluster = config['queue_small']+" -o "+basedir+"/LOG/index_dp_vcf_{chrm}.out -e "+basedir+"/LOG/index_dp_vcf_{chrm}.err -N index_dp_vcf_{chrm}",
        vcf_file = data_dir+'/VCF_FILTERED/'+prefix+"{chrm}"+suffix+"_dp6.vcf"
    run:
        shell("gunzip -f {input.vcf_file}")
        shell("bgzip -f {params.vcf_file}")
        shell("bcftools index --tbi {input.vcf_file}")

rule add_ANC_to_VCF:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF_FILTERED/baylor_phased_chr{chrm}_dp6.vcf.gz"
    output:
        vcf_file = data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc.vcf.gz"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/add_ANC_to_VCF_{chrm}.out -N add_ANC_to_VCF_{chrm}",
        vcf_file = data_dir+"/VCF_FILTERED/baylor_phased_chr{chrm}_dp6.vcf",
        vcf_file_out = data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc.vcf"
    run:
        shell("gunzip -c {input.vcf_file} > {params.vcf_file}")
        shell(homedir+"/scripts/add-ANC-to-vcf_new.py -g --in {params.vcf_file} --out {params.vcf_file_out} --genomedata "+genomedata_path)
        shell("bgzip -f {params.vcf_file_out}")

rule filter_ANC_only_VCF:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc.vcf.gz"
    output:
        vcf_file = data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc_f.vcf.gz"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/filter_ANC_{chrm}.out -N filter_ANC_{chrm}"
    run:
        shell(VENV_DIR+"/bcftools view -i 'AA!=\".\" & AA!=\"-\" & AA!=\"N\"' {input.vcf_file} | \
            "+VENV_DIR+"/bgzip -c > {output.vcf_file}")

rule download_database_snpEff:
    output:
        db = snpEff_database+"/"+snpEff_human_db+".txt"
    params:
        cluster = config['queue_medium']+" -o  "+basedir+"/LOG/snpEff_download_db.out -j oe -N snpEff_download"
    run:
        shell(snpEff+" download "+snpEff_human_db+" -verbose \
            -dataDir "+snpEff_database)
        shell("touch {output.db}")

# rule download_dbNSFT_snpSift:
#     output:
#         db = snpEff_database+"/dbNSFT.txt"
#     params:
#         cluster = config['queue_medium']+" -o  "+basedir+"/LOG/snpEff_download_db.out -j oe -N snpEff_download"
#     run:
#         shell(snpSift+" -download "+snpEff_dbNSFT_db+" -verbose \
#             -dataDir "+snpEff_database)
#         shell("touch {output.db}")

rule download_dbsnp_snpEff:
    output:
        db = dbsnp_vcf
    params:
        cluster = config['queue_medium']+" -o  "+basedir+"/LOG/snpEff_download_dbsnp.out -j oe -N snpEff_download_dbsnp"
    run:
        shell("wget -O "+snpEff_database+"/dbsnp137.201206.vcf.gz "+snpEff_dbsnp_url)
        shell("gunzip "+snpEff_database+"/dbsnp137.201206.vcf.gz")

rule download_gwascatalog_snpEff:
    output:
        db = snpEff_database+"/gwascatalog.txt"
    params:
        cluster = config['queue_medium']+" -o  "+basedir+"/LOG/snpEff_download.out -j oe -N snpEff_download_gwascat"
    run:
        shell("wget -O "+snpEff_database+"/gwascatalog.txt "+snpEff_gwascatalog)

rule annotate_vcf_dbsnp_snpEff:
    input:
        vcf = data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc_f.vcf.gz",
        dbsnp = dbsnp_vcf
    output:
        vcf = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp.vcf.gz"
    params:
        cluster = config['queue_medium']+" -o  "+basedir+"/LOG/snpEff_id.{chrm}.out -j oe -N snpEff_id.{chrm}",
        vcf_in = data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_anc_f.vcf",
        vcf_out = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp.vcf"
    run:
        shell("gunzip -c {input.vcf} > {params.vcf_in}")
        shell(snpSift+ " \
            annotate \
            {input.dbsnp} \
            {params.vcf_in} > {params.vcf_out} -v ")
        shell("bgzip -f {params.vcf_out}")
        shell("rm -f {params.vcf_in}")

rule annotate_vcf_snpEff:
    input:
        snpEff_database+"/"+snpEff_human_db+".txt",
        vcf = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp.vcf.gz"
    output:
        vcf = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz"
    params:
        cluster = config['queue_medium']+" -o  "+basedir+"/LOG/snpEff.{chrm}.out -j oe -N snpEff.{chrm}",
        vcf = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp",
        vcf_ = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff"
    threads: 1
    run:
        shell("gunzip -c {params.vcf}.vcf.gz > {params.vcf}.vcf")
        shell(snpEff+" \
            -v "+ \
            snpEff_human_db+" \
            -stats {params.vcf_}.html \
            -csvStats {params.vcf_}.csv \
            -dataDir "+snpEff_database+" \
            {params.vcf}.vcf > {params.vcf_}.vcf")
        shell("bgzip -f {params.vcf_}.vcf")
        shell("rm -f {params.vcf}.vcf")

rule split_POP_samples:
    '''
    Split sample files by population
    '''
    input:
        sample_file = sample_file
    output:
        sample_file = basedir+"/samples/{POP}.sample",
    params:
        cluster = config['queue_small']+" -o "+basedir+"/LOG/sample_{POP}.out -j oe -N sample_{POP}"
    run:
        split_sample_list_per_pop(wildcards.POP, input_file=input.sample_file, output_file=output.sample_file)

rule index_vcf_ann:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz"
    output:
        vcf_file = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz.tbi"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/tbi_raw_{chrm}.out -N tbi_raw_{chrm}",
    run:
        shell("bcftools index --tbi {input.vcf_file}")

rule split_vcf_by_pop:
    """
    Split VCF by population using VCFTools
    """
    input:
        vcf = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz",
        tbi = data_dir+"/VCF_ANN/baylor_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz.tbi",
        pop = basedir+"/samples/{POP}.sample"
    output:
        vcf_pop = data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz"
    params:
        cluster = config['queue_smaller']+" -o "+basedir+"/LOG/split.vcf.{POP}.{chrm}.out -j oe -N split.{POP}.{chrm}"
    run:
        shell(VENV_DIR+"/vcftools \
            --gzvcf {input.vcf} \
            --keep {input.pop} \
            --recode --recode-INFO-all --mac 1 -c | "+ \
            VENV_DIR+"/bgzip \
            -c > {output.vcf_pop}")

rule split_vcf_by_pop_tbi:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz"
    output:
        vcf_file = data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz.tbi"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/tbi_{POP}_{chrm}.out -N tbi_{POP}_{chrm}",
    run:
        shell("bcftools index --tbi {input.vcf_file}")

rule concat_vcf_in_pop:
    """
    Concat VCF chrm in population using BCFTools
    """
    input:
        vcf = lambda wildcards: expand(data_dir+"/VCF_POP/"+wildcards.POP+"_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz", chrm=CHRS),
        tbi = lambda wildcards: expand(data_dir+"/VCF_POP/"+wildcards.POP+"_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz.tbi", chrm=CHRS),
    output:
        vcf = data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.vcf.gz"
    params:
        cluster = config['queue_smaller']+" -o "+basedir+"/LOG/concat.vcf.{POP}.out -j oe -N concat.{POP}"
    run:
        shell(VENV_DIR+"/bcftools concat {input.vcf} \
            -Oz -o {output.vcf}")

rule index_vcf_concat:
    '''
    '''
    input:
        vcf_file = data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.vcf.gz"
    output:
        vcf_file = data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.vcf.gz.tbi"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/tbi_{POP}.out -N tbi_{POP}",
    run:
        shell("bcftools index --tbi {input.vcf_file}")

rule daf_by_pop:
    """
    Compute DAF by population using VCFTools
    """
    input:
        vcf = data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.vcf.gz",
        tbi = data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.vcf.gz.tbi",
    output:
        frq = data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.daf.frq"
    params:
        cluster = config['queue_smaller']+" -o "+basedir+"/LOG/split.vcf.{POP}.out -j oe -N split.{POP}",
        frq = data_dir+"/VCF_POP/{POP}_phased_dp6_anc_f_dbsnp_snpeff.daf"
    run:
        shell(VENV_DIR+"/vcftools \
            --gzvcf {input.vcf} \
            --freq --derived \
            --out {params.frq}")

rule daf_by_pop_chr:
    """
    Compute DAF by population using VCFTools
    """
    input:
        vcf = data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz",
        tbi = data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz.tbi",
    output:
        frq = data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.daf.frq"
    params:
        cluster = config['queue_smaller']+" -o "+basedir+"/LOG/daf.{POP}.{chrm}.out -j oe -N daf.{POP}.{chrm}",
        frq = data_dir+"/VCF_POP/{POP}_phased_chr{chrm}_dp6_anc_f_dbsnp_snpeff.daf"
    run:
        shell(VENV_DIR+"/vcftools \
            --gzvcf {input.vcf} \
            --freq --derived \
            --out {params.frq}")

# rule update_anc_alg_1000g:
#     '''
#     '''
#     input:
#         fasta_bz2 = data_dir+"/reference/ancestral_alignments/human_ancestor_{chrm}.fa.bz2"
#     output:
#         fasta_gz = data_dir+"/reference/ancestral_alignments/human_ancestor_{chrm}.fa.gz"
#     params:
#         cluster = config['queue_big']+" -j oe -o "+basedir+"/LOG/update_anc_alg_{chrm}.out -N update_anc_alg_{chrm}",
#         fasta = data_dir+"/reference/ancestral_alignments/human_ancestor_{chrm}.fa"
#     run:
#         shell("bzcat {input.fasta_bz2} | \
#             sed 's,^>.*,>1,' | \
#             "+VENV_DIR+"/bgzip -c > {output.fasta_gz}")
#         shell(VENV_DIR+"/samtools faidx {output.fasta_gz}")

# rule fill_AA_to_VCF:
#     '''
#     '''
#     input:
#         vcf_file = data_dir+"/VCF_FILTERED/baylor_phased_chr{chrm}_dp6.vcf.gz",
#         ref_fasta_gz = data_dir+"/reference/ancestral_alignments/human_ancestor_{chrm}.fa.gz"
#     output:
#         vcf_file = data_dir+"/VCF_ANC/baylor_phased_chr{chrm}_dp6_AA.vcf.gz",
#     params:
#         cluster = config['queue_big']+" -j oe -o "+basedir+"/LOG/add_AA_{chrm}.out -N add_AA_{chrm}"
#     run:
#         shell("zcat {input.vcf_file} | " + \
#             VENV_DIR+"/fill-aa -a "+data_dir+"/reference/ancestral_alignments/human_ancestor_ | " + \
#             VENV_DIR+"/bgzip -c > {output.vcf_file}")

rule report_qc_DP:
    '''
    '''
    input:
        vcf_before = expand(data_dir+'/VCF/'+prefix+"{chrm}"+suffix+".vcf.gz", chrm=CHRS),
        tbi_passre = expand(data_dir+'/VCF/'+prefix+"{chrm}"+suffix+".vcf.gz.tbi", chrm=CHRS),
        vcf_after = expand(data_dir+'/VCF_FILTERED/'+prefix+"{chrm}"+suffix+"_dp6.vcf.gz", chrm=CHRS),
        tbi_vcf_after = expand(data_dir+'/VCF_FILTERED/'+prefix+"{chrm}"+suffix+"_dp6.vcf.gz.tbi", chrm=CHRS)
    output:
        dp_report = data_dir+'/'+prefix+"{chrm}"+suffix+"_dp6.csv"
    params:
        cluster = config['queue_small']+" -j oe -o "+basedir+"/LOG/dp_report.out -N dp_report"
    run:
        # get total variant count
        generate_dp_report(fil_bef=data_dir+'/VCF/'+prefix+"<CHRM>"+suffix+".vcf.gz", fil_after=data_dir+'/VCF_FILTERED/'+prefix+"<CHRM>"+suffix+"_dp6.vcf.gz", dp_report=dp_report, CHRS=CHRS)


#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/retired_reference/ancestral_alignments/README.ancestral_alignments## Get the files ready: compress by gzip and index by samtools faidx. Either repeat the
## following command for each file manually
#bzcat human_ancestor_1.fa.bz2 | sed 's,^>.*,>1,' | gzip -c > human_ancestor_1.fa.gz
#samtools faidx human_ancestor_1.fa.gz
#
## .. or use this loop (tested in bash shell)
#ls human_ancestor_*.fa.bz2 | while read IN; do
#OUT=`echo $IN | sed 's,bz2$,gz,'`
#CHR=`echo $IN | sed 's,human_ancestor_,, ; s,.fa.bz2,,'`
#bzcat $IN | sed "s,^>.*,>$CHR," | gzip -c > $OUT
#samtools faidx $OUT
#done
#
## After this has been done, the following command should return 'TACGTGGcTGCTCTCACACAT'
#samtools faidx human_ancestor_1.fa.gz 1:1000000-1000020
#
## Now the files are ready to use with fill-aa. Note that the VCF file
## should be sorted (see vcf-sort), otherwise the performance would be seriously
## affected.
#### zcat file.vcf | fill-aa -a human_ancestor_ 2>test.err | gzip -c >out.vcf.gz
