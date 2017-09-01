# This folder contains code created by Gerrit

1. Get multi-alleles - `qsub get_multi-alleles.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/multi_alleles_only` using the orignial unphased VCFs in `/spaces/gapw/diversity/gerrit/baylor_post_vqsr_clean`
2. Remove INDELs from multi-allele sites. `qsub remove_indels.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/no_indels` using the multi-allele only files in `/spaces/gapw/diversity/gerrit/multi_alleles_only`
3. Getting MAFs from TrypanoGEN unphased and phased VCFs. `qsub get_trypanogen_maf.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/trypanogen_maf/` using the multi-allele only files in `/dataB/popdata/gapw/GAPW_DATA/trypanogen/UNPHASED/*.trypanogen_post_vqsr_clean.norm.vcf.gz` and `/dataB/popdata/gapw/GAPW_DATA/trypanogen/Eagle.trypanogen.*.vcf.gz

## Some additional one liners
*Create one column file for IBT to be filltered*
* `cat /spaces/gapw/diversity/filter/related.remove | awk '{print "^"$1}' > /spaces/gapw/diversity/filter/related.remove.mamana_ready`

*Get MAF annotation for SAHGP, AGVP and TrypanoGEN sets*
* `cat /spaces/gapw/diversity/dbs/agv3f.frq.frq | grep -v "CHR" | awk '{print $2"\t"$5}' | awk '{if($2!=0){print $0}}' > /spaces/gapw/diversity/dbs/agv3f.frq.frq.mamana_ready`
* `cat /spaces/gapw/diversity/dbs/sahgp_macs | awk '{print $1"\t"$3/30}' | sed "s/_/:/" > /spaces/gapw/diversity/dbs/sahgp_macs.mamana_ready`
* `grep -v "CHROM" /spaces/gapw/diversity/gerrit/trypanogen_maf/trypanogen.all.phased.frq  | awk '{print $1"_"$2"\t"$6}' | awk '{if($2!=0){print $0}}' | sed "s/.://" | sed "s/_/:/" > /spaces/gapw/diversity/dbs/trypanogen.all.phased.frq.mamana_ready`
