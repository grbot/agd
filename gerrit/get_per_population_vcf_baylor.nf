#!/usr/bin/env nextflow

in_dir = file(params.baylor_vcf_dir)
out = file(params.split_per_pop_vcfs_dir)

sample_to_population_mapping_dir = file(params.sample_to_population_mapping_dir)

out.mkdir()

pops = "BOT,BBC,BRN,BSZ,FNB,MAL,WGR".split(',')

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
	   file("${pop}.${chrom}.baylor.annotated.vcf.gz") into pop_chrom_file 

    """
    /opt/exp_soft/bcftools-1.6/install/bin/bcftools view -S ${sample_to_population_mapping_dir}/${pop}_sample_id_only.tsv ${in_dir}/Eagle.baylor.${chrom}.vcf.gz -O z -o ${pop}.${chrom}.baylor.vcf.gz
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

