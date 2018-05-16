## Multi-allele processing

1. Get multi-alleles - `qsub get_multi-alleles.qsub` - output files will be generated in 

`/space/gapw/diversity/gerrit/multi_alleles_only` using the orignial unphased VCFs in `/spaces/gapw/diversity/gerrit/baylor_post_vqsr_clean`

2. Remove INDELs from multi-allele sites. `qsub remove_indels.qsub` - output files will be generated in 

`/spaces/gapw/diversity/gerrit/no_indels` using the multi-allele only files in `/spaces/gapw/diversity/gerrit/multi_alleles_only`

3. Combine autosomes into one VCF, remove sites with high missingness, remove related individuals, remove sites if any sample has a read depth < 6.
```
module load bioinf
for i in {1..22}; do ls -1 /spaces/gapw/diversity/gerrit/no_indels/$i.baylor_post_vqsr_clean.multi_alleles_only.no_indels.vcf.gz; done > /spaces/gapw/diversity/gerrit/no_indels_combined/autosome.list 

bcftools concat -f /spaces/gapw/diversity/gerrit/no_indels_combined/autosome.list -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.vcf.gz -O z --threads 24

sed  "s/:/\t/" /spaces/gapw/diversity/filter/high_missing.snp > /spaces/gapw/diversity/filter/high_missing.for_bcftools.snp

bcftools view -T ^/spaces/gapw/diversity/filter/high_missing.for_bcftools.snp -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.vcf.gz -O z /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.vcf.gz

bcftools view --force-samples -S ^/spaces/gapw/diversity/filter/related.remove.mamana_ready -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.vcf.gz -O z /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.vcf.gz

bcftools view -i "MIN(FMT/DP)>=6" /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.vcf.gz -O z -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.vcf.gz

```
4. Creating some files to finally get "true multi-alleles"

  * Convert VCF to new format for better parsing (see `convert_vcf.py` for output format)
  
    ```
    ./convert_vcf.py -v /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.vcf.gz > /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.tsv
    ```

  * Get sites where REF = 0.000 but we have more than 3 alleles.
  
     ```
     cat /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.tsv | awk '{pos = index($3, ","); ac_ref=substr($3,1,pos-1);if(ac_ref == "0.000"){print $0}}'  | awk '{if(length($2) > 5){print $0}}' > /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000_more_than_tri-allele.tsv
     ```

  * Get sites where REF = 0.000 but we have  3 alleles.

    ```
    cat /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.tsv | awk '{pos = index($3, ","); ac_ref=substr($3,1,pos-1);if(ac_ref == "0.000"){print $0}}'  | awk '{if(length($2) == 5){print $0}}' > /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000_tri-allele.tsv
    ```
    
  * Get all REF = 0.000 sites for doing filtering with BCFtools using inclusion/exclusion (creating .bed file)
  
    ```
    cat /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000_tri-allele.tsv  /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000_more_than_tri-allele.tsv | awk '{split($1,coord,":"); print coord[1]"\t"(coord[2]-1)"\t"coord[2]}' > /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.bed
    ```
  
  * Get sites with "true" tri-allelles or more

    ```
    bcftools view -T ^/spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.bed  -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.vcf.gz -O z /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.vcf.gz
    ```


 * Get sites witn "non-true" ref=0.000 tri-alleles and some ref=0.000 quadri-alleles
 
    ```
    bcftools view -T /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.bed  -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.vcf.gz -O z /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.vcf.gz
    ```

5. Get per population VCFs.
```
for i in {"BBC","BOT","BRN","BSZ","FNB","MAL","WGR"}; do bcftools view -S /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i\_sample_id_only.tsv --min-ac=1 --force-samples  -O z -o /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.vcf.gz /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.vcf.gz; done
```

6. Get per population sample sample counts. Can be used to normalise per site using the table generated in 7.

```
for i in {"BBC","BOT","BRN","BSZ","FNB","MAL","WGR"}; do nr_samples=`zcat /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.vcf.gz |  head -n 200 | grep "^#CHROM" | cut -f10- | sed "s/\t/\n/g" | wc -l`; echo -n $i$'\t'; echo $nr_sites$'\t'$nr_samples$'\t'; done
```

7.Convert per population VCF to new format for better parsing.

```
for i in {"BBC","BOT","BRN","BSZ","FNB","MAL","WGR"}; do /home/gerrit/projects/agd/gerrit/convert_vcf.py -v /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.vcf.gz > /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.tsv ; done
```

8. Get table with allele count within populations. Hard coded `get_multi-allele_pop_stats.py` to point to `/spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.tsv`. Then ran:

```
./get_multi-allele_pop_stats.py
```

A table like this is generated:
```
# 1) Generates allele counts
# ---- Allele counts ----
# Population    Bi-allelic count        Tri-allelic count       Other-allelic count     Allele types
# 2) Generate allele scenario counts
# Population    Two common alleles and one rare allele  Two rare alleles and one common allele
```
To get the allele counts on the exome region change the hard coded path in `get_multi-allele_pop_stats.py` to point to `/spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exome.tsv` and run `./get_multi-allele_pop_stats.py` again.

9. Pull out exonic sites only

  * Downloaded ENSEMBL annotations from here  `http://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz`
    
  * Then pulled out exonic regions only

  ```
  grep "exon" Homo_sapiens.GRCh37.87.gtf | awk '{ if($1=='1'||$1=="2"||$1=="3"||$1=="4"||$1=="5"||$1=="6"||$1=="7"||$1=="8"||$1=="9"||$1=="10"||$1=="11"||$1=="12"||$1=="13"||$1=="14"||$1=="15"||$1=="16"||$1=="16"||$1=="18"||$1=="19"||$1=="20"||$1=="21"||$1=="22"||$1=="MT"||$1=="Y"||$1=="X"){print $0}}' > Homo_sapiens.GRCh37.87.exon.gtf
  ```

  * Converted to .bed
     
    ```
    ~/software/bedops/bin/gtf2bed < Homo_sapiens.GRCh37.87.exon.gtf > Homo_sapiens.GRCh37.87.exon.bed
    ```
    
  * Copied to AGD folder so that everyone can access it
 
    ```
    cp /global/chpdes/gerrit/prep_exonic_non_exonic_bead_pools/Homo_sapiens.GRCh37.87.exon.bed    /spaces/gapw/diversity/filter/
    ```
  
  * Get sites with "true" tri-allelles or more, exonic regions only. Also convert to tab format with `convert_vcf.py`.

    ```
    bcftools view -T /spaces/gapw/diversity/filter/Homo_sapiens.GRCh37.87.exon.bed -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exon.vcf.gz -O z /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.vcf.gz
    /home/gerrit/projects/agd/gerrit/convert_vcf.py -v /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exon.vcf.gz > /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exon.tsv
    ```

  * Get sites witn "non-true" ref=0.000 tri-alleles and some ref=0.000 quadri-alleles, exonic regions only. Also convert to tab format with `convert_vcf.py`.
  
    ```
    bcftools view -T /spaces/gapw/diversity/filter/Homo_sapiens.GRCh37.87.exon.bed -o /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.exon.vcf.gz -O z /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.vcf.gz
    /home/gerrit/projects/agd/gerrit/convert_vcf.py -v /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.exon.vcf.gz > /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.ref_0.000.exon.tsv
    ```

  * Get per population VCFs, exonic regions only. Also convert to tab format with `convert_vcf.py`.
  
    ```
    for i in {"BBC","BOT","BRN","BSZ","FNB","MAL","WGR"}; do bcftools view -S /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i\_sample_id_only.tsv --min-ac=1 --force-samples  -O z -o /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exon.vcf.gz /spaces/gapw/diversity/gerrit/no_indels_combined/all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exon.vcf.gz; done
    for i in {"BBC","BOT","BRN","BSZ","FNB","MAL","WGR"}; do /home/gerrit/projects/agd/gerrit/convert_vcf.py -v /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exon.vcf.gz > /spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs/$i.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.exon.tsv ; do
    ``` 

## MAF processing

1. Getting MAFs from TrypanoGEN unphased and phased VCFs. `qsub get_trypanogen_maf.qsub` - output files will be generated in `/spaces/gapw/diversity/gerrit/trypanogen_maf/` using the multi-allele only files in `/dataB/popdata/gapw/GAPW_DATA/trypanogen/UNPHASED/*.trypanogen_post_vqsr_clean.norm.vcf.gz` and `/dataB/popdata/gapw/GAPW_DATA/trypanogen/Eagle.trypanogen.*.vcf.gz`

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
 
