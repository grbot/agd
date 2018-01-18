## Multi-allele processing

1. Get multi-alleles - `qsub get_multi-alleles.qsub` - output files will be generated in `/space/gapw/diversity/gerrit/multi_alleles_only` using the orignial unphased VCFs in `/spaces/gapw/diversity/gerrit/baylor_post_vqsr_clean`
2. Remove INDELs from multi-allele sites. `qsub remove_indels.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/no_indels` using the multi-allele only files in `/spaces/gapw/diversity/gerrit/multi_alleles_only`


## MAF processing

1. Getting MAFs from TrypanoGEN unphased and phased VCFs. `qsub get_trypanogen_maf.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/trypanogen_maf/` using the multi-allele only files in `/dataB/popdata/gapw/GAPW_DATA/trypanogen/UNPHASED/*.trypanogen_post_vqsr_clean.norm.vcf.gz` and `/dataB/popdata/gapw/GAPW_DATA/trypanogen/Eagle.trypanogen.*.vcf.gz

## Some additional one liners
*Create one column file for IBT to be filltered*
* `cat /spaces/gapw/diversity/filter/related.remove | awk '{print "^"$1}' > /spaces/gapw/diversity/filter/related.remove.mamana_ready`

*Get MAF annotation for SAHGP, AGVP and TrypanoGEN sets*
* `cat /spaces/gapw/diversity/dbs/agv3f.frq.frq | grep -v "CHR" | awk '{print $2"\t"$5}' | awk '{if($2!=0){print $0}}' > /spaces/gapw/diversity/dbs/agv3f.frq.frq.mamana_ready`
* `cat /spaces/gapw/diversity/dbs/sahgp_macs | awk '{print $1"\t"$3/30}' | sed "s/_/:/" > /spaces/gapw/diversity/dbs/sahgp_macs.mamana_ready`
* `grep -v "CHROM" /spaces/gapw/diversity/gerrit/trypanogen_maf/trypanogen.all.phased.frq  | awk '{print $1"_"$2"\t"$6}' | awk '{if($2!=0){print $0}}' | sed "s/.://" | sed "s/_/:/" > /spaces/gapw/diversity/dbs/trypanogen.all.phased.frq.mamana_ready`

## Some QC checks

1. Get sites that have `ALT=.` and also count the total number of sites in the merged phased set -`nextflow -log nextflow.log run -w /spaces/gapw/diversity/gerrit/nextflow/workdir -c /home/gerrit/projects/agd/gerrit/nextflow.config /home/gerrit/projects/agd/gerrit/get_unidentified_alt_from_merged_phased.nf -profile pbs`


## Calculating per population SNP stats
1. Test nextflow script - `nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/gapw/nextflow/workdir -c /home/gerrit/code/agd/gerrit/nextflow.config /home/gerrit/code/agd/gerrit/test.nf -profile pbs` 

2. Get per population VCFs - `nextflow -log nextflow.log run -w /spaces/gapw/diversity/gerrit/nextflow/workdir -c /home/gerrit/projects/agd/gerrit/nextflow.config /home/gerrit/projects/agd/gerrit/get_per_population_vcf.nf -profile pbs` (This ran successfully at Wits).

3. Concatenate per population chromosome VCFs - `nextflow -log nextflow.log run -w /spaces/gapw/diversity/gerrit/nextflow/workdir -c /home/gerrit/projects/agd/gerrit/nextflow.config /home/gerrit/projects/agd/gerrit/concat_per_population_vcf.nf -profile pbs`

4. Annotating with VEP. This was done on CBIO because the VEP container and db setup is there.

(process 10 populations)

```
qsub -I -q dev -l nodes=1:ppn=64 -l walltime=168:00:00

for i in {"BBC","BOT","BRN","BSZ","CIV","DRC","FNB","MAL","SOT","SSG"}; do echo $i; done | /home/gerrit/soft/parallel-20171022/install/bin/parallel singularity exec -H /global5/scratch/gerrit:/scratch /global5/scratch/gerrit/singularity-containers/ensemblorg_ensembl-vep.img perl /home/vep/src/ensembl-vep/vep --offline --cache --dir /scratch/dbs/vep/ --species homo_sapiens --assembly GRCh37 --format vcf --af_gnomad --tab -i /scratch/projects/gapw/concat_per_pop_vcfs/{}/{}.Eagle.merged.all_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.vcf.gz -o /scratch/projects/gapw/per_pop_veps/{}.Eagle.merged.all_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.tsv --force_overwrite --sift b --polyphen b --humdiv --regulatory --hgvs --symbol --uniprot --af --max_af --af_1kg -af_gnomad --fork 6
```
(process 4 populations)

```
qsub -I -l nodes=1:ppn=30 -l walltime=168:00:00

for i in {"UBS","UNS","WGR","XHS"}; do echo $i; done | /home/gerrit/soft/parallel-20171022/install/bin/parallel singularity exec -H /global5/scratch/gerrit:/scratch /global5/scratch/gerrit/singularity-containers/ensemblorg_ensembl-vep.img perl /home/vep/src/ensembl-vep/vep --offline --cache --dir /scratch/dbs/vep/ --species homo_sapiens --assembly GRCh37 --format vcf --af_gnomad --tab -i /scratch/projects/gapw/concat_per_pop_vcfs/{}/{}.Eagle.merged.all_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.vcf.gz -o /scratch/projects/gapw/per_pop_veps/{}.Eagle.merged.all_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.tsv --force_overwrite --sift b --polyphen b --humdiv --regulatory --hgvs --symbol --uniprot --af --max_af --af_1kg -af_gnomad --fork 6
```
 
