# This folder contains code created by Gerrit

1. Get multi-alleles - `qsub get_multi-alleles.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/multi_alleles_only` using the orignial unphased VCFs in `/spaces/gapw/diversity/gerrit/baylor_post_vqsr_clean`
2. Remove INDELs from multi-allele sites. `qsub remove_indels.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/no_indels` using the multi-allele only files in `/spaces/gapw/diversity/gerrit/multi_alleles_only`
