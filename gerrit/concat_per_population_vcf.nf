#!/usr/bin/env nextflow

in_dir = file(params.split_per_pop_vcfs_dir)
out = file(params.concat_per_pop_vcfs_dir)

pops = "BBC,BOT,BRN,BSZ,CIV,DRC,FNB,MAL,SOT,SSG,UBS,UNS,WGR,XHS".split(',')

pops_ch = Channel.from(pops)

process concatPerPopulationVCF {
     tag { "${params.project_name}.${pop}.cPPV" }
     
     publishDir "${out}/${pop}", mode: 'copy', overwrite: false

     input:
 	 val(pop) from pops_ch 

    output:
	   file("${pop}.Eagle.merged.all_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.vcf.gz") into pop_concat_file 

    """
    for i in {1..22}; do ls -1 ${in_dir}/${pop}/*.\$i\\_*.vcf.gz; done > vcf.list
    /opt/exp_soft/bioinf/bin/bcftools concat -f vcf.list -O z -o ${pop}.Eagle.merged.all_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_noALT.vcf.gz
    """

}
pop_concat_file.subscribe { println it }

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

