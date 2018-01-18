# This folder contains code created by Gerrit

1. Get multi-alleles - `qsub get_multi-alleles.qsub` - output files will be generated in `/space/gapw/diversity/gerrit/multi_alleles_only` using the orignial unphased VCFs in `/spaces/gapw/diversity/gerrit/baylor_post_vqsr_clean`
2. Remove INDELs from multi-allele sites. `qsub remove_indels.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/no_indels` using the multi-allele only files in `/spaces/gapw/diversity/gerrit/multi_alleles_only`
3. Getting MAFs from TrypanoGEN unphased and phased VCFs. `qsub get_trypanogen_maf.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/trypanogen_maf/` using the multi-allele only files in `/dataB/popdata/gapw/GAPW_DATA/trypanogen/UNPHASED/*.trypanogen_post_vqsr_clean.norm.vcf.gz` and `/dataB/popdata/gapw/GAPW_DATA/trypanogen/Eagle.trypanogen.*.vcf.gz

## Some additional one liners
*Create one column file for IBT to be filltered*
* `cat /spaces/gapw/diversity/filter/related.remove | awk '{print "^"$1}' > /spaces/gapw/diversity/filter/related.remove.mamana_ready`

*Get MAF annotation for SAHGP, AGVP and TrypanoGEN sets*
* `cat /spaces/gapw/diversity/dbs/agv3f.frq.frq | grep -v "CHR" | awk '{print $2"\t"$5}' | awk '{if($2!=0){print $0}}' > /spaces/gapw/diversity/dbs/agv3f.frq.frq.mamana_ready`
* `cat /spaces/gapw/diversity/dbs/sahgp_macs | awk '{print $1"\t"$3/30}' | sed "s/_/:/" > /spaces/gapw/diversity/dbs/sahgp_macs.mamana_ready`
* `grep -v "CHROM" /spaces/gapw/diversity/gerrit/trypanogen_maf/trypanogen.all.phased.frq  | awk '{print $1"_"$2"\t"$6}' | awk '{if($2!=0){print $0}}' | sed "s/.://" | sed "s/_/:/" > /spaces/gapw/diversity/dbs/trypanogen.all.phased.frq.mamana_ready`

# Some checks.

1. Get sites that have `ALT=.` and also count the total number of sites in the merged phased set -`nextflow -log nextflow.log run -w /spaces/gapw/diversity/gerrit/nextflow/workdir -c /home/gerrit/projects/agd/gerrit/nextflow.config /home/gerrit/projects/agd/gerrit/get_unidentified_alt_from_merged_phased.nf -profile pbs`


# Calculating per population SNP stats
1. Test nextflow script - `nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/gapw/nextflow/workdir -c /home/gerrit/code/agd/gerrit/nextflow.config /home/gerrit/code/agd/gerrit/test.nf -profile pbs` 

2. Get per population VCFs - `nextflow -log nextflow.log run -w /spaces/gapw/diversity/gerrit/nextflow/workdir -c /home/gerrit/projects/agd/gerrit/nextflow.config /home/gerrit/projects/agd/gerrit/get_per_population_vcf.nf -profile pbs` (This ran successfully at Wits).

3. Concatenate per population chromosome VCFs - `nextflow -log nextflow.log run -w /spaces/gapw/diversity/gerrit/nextflow/workdir -c /home/gerrit/projects/agd/gerrit/nextflow.config /home/gerrit/projects/agd/gerrit/concat_per_population_vcf.nf -profile pbs`
 
