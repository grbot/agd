# This folder contains code created by Gerrit

1. Get multi-alleles - `qsub get_multi-alleles.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/multi_alleles_only` using the orignial unphased VCFs in `/spaces/gapw/diversity/gerrit/baylor_post_vqsr_clean`
2. Remove INDELs from multi-allele sites. `qsub remove_indels.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/no_indels` using the multi-allele only files in `/spaces/gapw/diversity/gerrit/multi_alleles_only`
3. Getting MAFs from TrypanoGEN unphased and phased VCFs.. `qsub get_trypanogen_maf.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/trypanogen_maf/` using the multi-allele only files in `/dataB/popdata/gapw/GAPW_DATA/trypanogen/UNPHASED/*.trypanogen_post_vqsr_clean.norm.vcf.gz` and `/dataB/popdata/gapw/GAPW_DATA/trypanogen/Eagle.trypanogen.*.vcf.gz`
