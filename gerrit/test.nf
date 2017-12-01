#!/usr/bin/env nextflow
out_path = file(params.out)

out_path.mkdir()

chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y".split(',')

process loopChromosomes {
     tag { "${params.project_name}.${chr}.lC" }
     
     publishDir "${out_path}/${chr}", mode: 'copy', overwrite: false

     input:
	 each chr from chromosomes  

    output:
	   file("chr.txt") into chr_file 

    """
    echo "${chr}" > chr.txt
    """
}

chr_file.subscribe { println it }

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

