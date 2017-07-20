#!/usr/bin/env nextflow

/*
 * Authors:
 *      Mamana Mbiyavanga
 *
 *  On behalf of the H3ABionet Consortium
 *  2017
 *
 *
 * Description  : Nextflow pipeline for ...
 *
*/

//---- General definitions --------------------------------------------------//

CHRMS = params.chromosomes.split(',')

// All POP
def POPS_ALL = []
params.POPS.each { entry->
    POPS_ALL.addAll(entry.value.split(','))
}
println "Project : $workflow.projectDir"
//println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Populations in use: ${POPS_ALL.join(", ")}"
println "Chromosomes used: ${CHRMS.join(', ')}"



//// Help functions

def split_sample_list_per_pop(POP_, input_file, output_file){
    '''
    Read integrated_call_samples_v3.20130502.ALL.panel and split it by population
    '''
    data = []
    println("Extracting "+POP_+' to '+output_file)
    new File(input_file).eachLine { line, nline ->
        if(nline > 1){
            line = line.trim().split()
            POP = line[1].trim()
            if(POP == POP_){
                data.add(line.join('\t')+'\n')
            }
        }
    }
    data = data.join(' ')
    myFile = file(output_file)
    myFile.text = data
    return data
}

// Create a channel for initial data which is in chromosomes
datas = []
CHRMS.each { chromosome ->
    datas << [chromosome, file(sprintf(params.data, chromosome))]
}
vcf_data = Channel.from(datas)


'''
Step 1: Filter VCF by coverage depth using bcftools
'''
process filter_vcf_by_depth {
    tag "dp6_${chrm}"
    memory { 4.GB * task.attempt }
    time = { 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_FILTERED/", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_FILTERED/", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file) from vcf_data
    output:
        set val(chrm), file("${params.prefix_new}${chrm}_dp6.vcf.gz"), file("${params.prefix_new}${chrm}_dp6.vcf.gz.tbi") into vcf_6depth
    script:
        """
        bcftools view -i 'DP>6' ${vcf_file} | \
            bgzip -c > ${params.prefix_new}${chrm}_dp6.vcf.gz
        bcftools index --tbi -f ${params.prefix_new}${chrm}_dp6.vcf.gz
        """
}


'''
Step 2: Annotate VCF for Ancestral Alle (AA) using in-house python script
'''
vcf_6depth.into { vcf_6depth; vcf_6depth__anc } // Duplicate channel so that it can be used multiple times
process add_ANC_to_VCF {
    echo true
    tag "AA_${chrm}"
    memory { 8.GB * task.attempt }
    time = { 4.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_file), file(vcf_file_tbi) from vcf_6depth__anc
    output:
        set val(chrm), file("${vcf_file_out}.gz"), file("${vcf_file_out}.gz.tbi") into vcf_anc
    script:
        vcf_file_out = "${params.prefix_new}${chrm}_dp6_anc.vcf"
        """
        gunzip -c ${vcf_file} > ${vcf_file.baseName}
        ${params.python_27_env} ${params.homedir}/scripts/add-ANC-to-vcf_new.py -g --in ${vcf_file.baseName} --out ${vcf_file_out} --genomedata ${params.genomedata_path}
        bgzip -f ${vcf_file_out}
        bcftools index --tbi -f ${vcf_file_out}.gz
        rm -f ${vcf_file.baseName}
        """
}


'''
Step 3: Filter sites with AA
'''
vcf_anc.into { vcf_anc; vcf_anc__only}
process add_ANC_to_VCF_only {
    echo true
    tag "AA_only_${chrm}"
    memory { 2.GB * task.attempt }
    time = { 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_ANC/", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_dp6_anc), file(vcf_dp6_anc_tbi) from vcf_anc__only
    output:
        set val(chrm), file("${params.prefix_new}${chrm}_dp6_anc_f.vcf.gz"), file("${params.prefix_new}${chrm}_dp6_anc_f.vcf.gz.tbi") into vcf_anc_f
    script:
        vcf_dp6_anc_out = "${params.prefix_new}${chrm}_dp6_anc_f.vcf.gz"
        """
        bcftools view -i 'AA!="." & AA!="-" & AA!="N"' ${vcf_dp6_anc} | bgzip -c > ${vcf_dp6_anc_out}
        bcftools index --tbi -f ${vcf_dp6_anc_out}
        """
}

// TODO Download dbSNP database if not exists
//vcf_anc_f.into { vcf_anc_f; vcf_anc_f__dbsnp}
//process download_snpeff_db {
//    echo true
//    tag "dbSNP_${chrm}"
//    memory { 8.GB * task.attempt }
//    time = { 6.hour * task.attempt }
//    publishDir "${params.work_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
//    input:
//        set val(chrm), file(vcf_dp6_anc_f) from vcf_anc_f__dbsnp
//    output:
//        set val(chrm), file("${vcf_out}.gz") into annotate_dbsnp_snpeff_all
//    script:
//        """
//        download
//        index
//
//        """
//}


'''
Step 4:
'''
vcf_anc_f.into { vcf_anc_f; vcf_anc_f__dbsnp}
process annotate_dbsnp_snpeff {
    echo true
    tag "dbSNP_${chrm}"
    memory { 8.GB * task.attempt }
    time = { 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_dp6_anc_f) from vcf_anc_f__dbsnp
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_dbsnp_snpeff_all
    script:
        vcf_out = "${params.prefix_new}${chrm}_dp6_anc_f_dbsnp.vcf"
        """
        SnpSift \
            annotate \
            ${params.dbsnp_vcf} \
            ${vcf_dp6_anc_f} > ${vcf_out} -v
        bgzip -f ${vcf_out}
        """
}


'''
Step 5
'''
process annotate_snpeff {
    echo true
    tag { "snpEff_${chrm}" }
    memory { 8.GB * task.attempt }
    time = { 6.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_ANN/", overwrite: true, mode:'symlink'
    input:
        set val(chrm), file(vcf_dp6_anc_f_dbsnp) from annotate_dbsnp_snpeff_all
    output:
        set val(chrm), file("${vcf_out}.gz") into annotate_snpeff_all
    script:
        vcf_in = "${params.prefix_new}${chrm}_dp6_anc_f_dbsnp_snpeff"
        vcf_out = "${params.prefix_new}${chrm}_dp6_anc_f_dbsnp_snpeff.vcf"
        """
        ${params.snpEff} \
            ${params.snpEff_human_db} \
            -stats ${vcf_in}.html \
            -csvStats ${vcf_in}.csv \
            -dataDir ${params.snpEff_database} \
            ${vcf_dp6_anc_f_dbsnp} > ${vcf_out} -v -nodownload
        bgzip -f ${vcf_out}
        """
}


'''
Step 6: Generate
'''
process split_POP_samples {
    echo true
    tag { "split_POP_samples_${POP}" }
    memory { 2.GB  * task.attempt }
    time = { 1.hour * task.attempt }
    publishDir "${params.work_dir}/samples/", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/samples/", overwrite: true, mode:'symlink'
    input:
        val POP from POPS_ALL
    output:
        set val(POP), file("${POP}.sample") into split_POP_samples_all
    script:
        """
        grep ${POP} ${params.sample_file} | cut -f1- > ${POP}.sample
        """
}

'''
Step 7: Split vcf per population
'''
annotate_snpeff_all.into { annotate_snpeff_all; annotate_snpeff_all__split_pop}
pop_datas = []
annotate_snpeff_all__split_pop_list = annotate_snpeff_all__split_pop.toSortedList().val
split_POP_samples_all.into {split_POP_samples_all; split_POP_samples_all__split_pop}
split_POP_samples_all__split_pop.toSortedList().val.each { POP, sample_file ->
    annotate_snpeff_all__split_pop_list.each { chrm, vcf_file ->
        pop_datas << [POP, sample_file, chrm, vcf_file]
    }
}
annotate_snpeff_all__split_pop_cha = Channel.from(pop_datas)

process split_vcf_per_pop {
    echo true
    tag { "split_vcf_${POP}_${chrm}" }
    memory { 2.GB * task.attempt }
    time = { 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/${POP}", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_POP/${POP}", overwrite: true, mode:'symlink'
    input:
        set val(POP), file(POP_sample_file), val(chrm), file(chrm_vcf_file) from annotate_snpeff_all__split_pop_cha
    output:
        set val(POP), val(chrm), file(vcf_out) into split_vcf_per_pop
    script:
        vcf_out = "${POP}${chrm}_dp6_anc_f_dbsnp_snpeff.vcf.gz"
        """
        vcftools \
            --gzvcf ${chrm_vcf_file} \
            --keep ${POP_sample_file} \
            --recode --recode-INFO-all -c | \
        bgzip -c > ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}

'''
Step: DAF per population chromosomes
'''
split_vcf_per_pop.into { split_vcf_per_pop; split_vcf_per_pop__daf_by_chrm_pop }
process daf_by_chrm_pop {
    echo true
    tag { "daf_${POP}_${chrm}" }
    memory { 2.GB * task.attempt }
    time = { 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/DAF", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_POP/DAF", overwrite: true, mode:'symlink'
    input:
        set val(POP), val(chrm), file(vcf_pop_chrm) from split_vcf_per_pop__daf_by_chrm_pop
    output:
        set val(POP), val(chrm), file("${vcf_out}.frq") into daf_by_chrm_pop
    script:
        vcf_out = "${POP}_${chrm}_dp6_anc_f_dbsnp_snpeff.daf"
        """
        vcftools \
            --gzvcf ${vcf_pop_chrm} \
            --freq --derived \
            --out ${vcf_out}
        """
}
daf_by_chrm_pop.into{ daf_by_chrm_pop; daf_by_chrm_pop_sub}
daf_by_chrm_pop_sub.subscribe {
    println "|-- Finished for ${it[-1][-1]}"
}

'''
Step: Concatenate chromosome VCFs into one  
'''
split_vcf_per_pop.into { split_vcf_per_pop; split_vcf_per_pop__merge_vcf_pop}
merge_vcf_pop_list = [:]
split_vcf_per_pop__merge_vcf_pop.toSortedList().val.each { POP, chrm, vcf_chrm_pop ->
    if ( !(POP in merge_vcf_pop_list.keySet()) ) {
        merge_vcf_pop_list[POP] = [POP]
    }
    else{
        merge_vcf_pop_list[POP] << vcf_chrm_pop
    }
}
merge_vcf_pop_cha = Channel.from(merge_vcf_pop_list.values())
merge_vcf_pop_cha.into { merge_vcf_pop_cha; merge_vcf_pop__merge}
process merge_vcf_pop {
    echo true
    tag { "merge_vcf_pop_${POP}" }
    memory { 5.GB * task.attempt }
    time = { 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/${POP}", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_POP/${POP}", overwrite: true, mode:'symlink'
    input:
        val datas from merge_vcf_pop__merge
    output:
        set val(POP), file(vcf_out) into merge_vcf_pop
    script:
        POP = datas[0]
        vcf_data = datas[1..-1].join(' ')
        vcf_out = "${POP}_dp6_anc_f_dbsnp_snpeff.vcf.gz"
        """
        bcftools concat ${vcf_data} -Oz -o ${vcf_out}
        bcftools index --tbi -f ${vcf_out}
        """
}


merge_vcf_pop.into { merge_vcf_pop; merge_vcf_pop__daf}
process daf_by_pop {
    echo true
    tag { "daf_${POP}" }
    memory { 5.GB * task.attempt }
    time = { 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/DAF", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_POP/DAF", overwrite: true, mode:'symlink'
    input:
        set val(POP), file(vcf_pop) from merge_vcf_pop__daf
    output:
        set val(POP), file("${vcf_out}.frq") into daf_by_pop
    script:
        vcf_out = "${POP}_dp6_anc_f_dbsnp_snpeff.daf"
        """
        vcftools \
            --gzvcf ${vcf_pop} \
            --freq --derived \
            --out ${vcf_out}
        """
}


'''
Extract number of singletons per population
'''
merge_vcf_pop.into { merge_vcf_pop; merge_vcf_pop__singl}
process singl_by_pop {
    echo true
    tag { "singl_${POP}" }
    memory { 5.GB * task.attempt }
    time = { 2.hour * task.attempt }
    publishDir "${params.work_dir}/VCF_POP/SINGL", overwrite: true, mode:'symlink'
    publishDir "${params.group_dir}/VCF_POP/SINGL", overwrite: true, mode:'symlink'
    input:
        set val(POP), file(vcf_pop) from merge_vcf_pop__singl
    output:
        set val(POP), file("${POP}.singletons.per.sample") into singl_by_pop
    script:
        vcf_out = "${POP}_dp6_anc_f_dbsnp_snpeff_derived"
        """
        vcftools \
            --gzvcf ${vcf_pop} \
            --singletons --derived \
            --out ${vcf_out}
        grep ${POP} ${params.sample_file} | cut -f1 > ${POP}.sample
        echo "singletons" > ${POP}.singleton.count
        while read SAMPLE; do 
            grep \$SAMPLE ${vcf_out}.singletons | wc -l;
        done < ${POP}.sample >> ${POP}.singleton.count
        echo "sample" > ${POP}_.sample | cat ${POP}.sample >> ${POP}_.sample
        paste ${POP}_.sample ${POP}.singleton.count > ${POP}.singletons.per.sample 
        """
}


// To run different steps
// if(params.target=='vcf_6depth'){
//     target = vcf_6depth
// }
// else if(params.target=='vcf_6depth_tbi'){
//     target = vcf_6depth_tbi
// }


workflow.onComplete {
    def subject = 'My pipeline execution'
    def recipient = 'mypandos@gmail.com'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

//annotate_snpeff_all.subscribe {
//    println "|-- Finished for ${it[1]}"
//}
