#!/usr/bin/env nextflow

in_dir = file(params.merged_annotated_vcf_dir)
out = file(params.split_per_pop_vcfs_dir)

sample_to_population_mapping_dir = file(params.sample_to_population_mapping_dir)

out.mkdir()

pops = "BBC,BOT,BRN,BSZ,CIV,DRC,FNB,MAL,SOT,SSG,UBS,UNS,WGR,XHS".split(',')

chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')

pops_ch = Channel.from(pops)

pops_ch.into { pops_ch_p1; pops_ch_p2 }

process loopPopulations {
     tag { "${params.project_name}.${pop}.lP" }
     
     publishDir "${out}/${pop}", mode: 'copy', overwrite: false

     input:
	 val(pop) from pops_ch_p1  

    output:
	   file("${pop}.txt") into pop_file 

    """
    echo "${pop}" > ${pop}.txt
    """
}

process getPerPopulationVCF {
     tag { "${params.project_name}.${pop}.${chrom}.gPPV" }
     
     publishDir "${out}/${pop}", mode: 'copy', overwrite: false

     input:
 	 val(pop) from pops_ch_p2 
	 each chrom from chroms  

    output:
	   file("${pop}.Eagle.merged.${chrom}_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.vcf.gz") into pop_chrom_file 

    """
    /opt/exp_soft/bioinf/bin/bcftools view -S ${sample_to_population_mapping_dir}/${pop}_sample_id_only.tsv ${in_dir}//Eagle.merged.${chrom}_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.vcf.gz -O z -o ${pop}.Eagle.merged.${chrom}_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.vcf.gz
    """

}
pop_chrom_file.subscribe { println it }

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

