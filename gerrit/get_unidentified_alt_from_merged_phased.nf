#!/usr/bin/env nextflow

in_dir = file(params.merged_phased_annotated_dir)
out = file(params.post_process_merged_phased_annotated)

out.mkdir()

chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')

process getUnIdentifiedAlt {
     tag { "${params.project_name}.${chrom}.gUIA" }
     
     publishDir "${out}/", mode: 'copy', overwrite: false

     input:
	 each chrom from chroms  

    output:
	   file("${chrom}.merged.phased.annotated.unidentified_alt.sites") into unidentified_alt_file 

    """
    zcat ${in_dir}/Eagle.merged.${chrom}_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss.vcf.gz | grep -v "^#" | awk '{if(\$5=="."){print \$1":"\$2}}' > ${chrom}.merged.phased.annotated.unidentified_alt.sites 
    """

}

process getTotalNrSites {
     tag { "${params.project_name}.${chrom}.gTNS" }
     
     publishDir "${out}/", mode: 'copy', overwrite: false

     input:
	 each chrom from chroms  

    output:
	   file("${chrom}.merged.phased.annotated.site_count.txt") into site_count_file 

    """
    zcat ${in_dir}/Eagle.merged.${chrom}_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss.vcf.gz | grep -v "^#" | wc -l > ${chrom}.merged.phased.annotated.site_count.txt 
    """
}

unidentified_alt_file.subscribe { println it }
site_count_file.subscribe { println it }

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

