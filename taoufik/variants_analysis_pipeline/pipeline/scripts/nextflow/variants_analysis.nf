#!/usr/bin/env nextflow
//params.afmin = 0.0
//params.gnomad = 0.0
in_dir = file(params.datavcf)
in_dir_freq = file(params.datafreq)
out_dir_stats = file("${params.output_dir}/NovelStats")
out_dir_stats.mkdirs()
out_dir = file("${params.output_dir}/NovelVariants")
out_dir.mkdirs()
out_dir_freq = file("${params.output_dir}/VariantsFreq")
out_dir_freq.mkdirs()
out_dir_impact = file("${params.output_dir}/VariantsImpact")
out_dir_impact.mkdirs()
//out_dir_sdt = file("${params.homeagd}/sdt")
//out_dir_sdt.mkdirs()
CHRMS = params.chromosomes.split(',')
POPS= params.pops.split(',')
//chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')
pops_ch = Channel.from(POPS)
vcfs = new ArrayList<String>();
ffreqcounts = new ArrayList<String>();
sdtouts = new ArrayList<String>()
Map<String,List> mappopvcfs  = new HashMap<String,List>();
Map<String,List> mappopffreqcounts  = new HashMap<String,List>();
println """\
             NOVEL VARIANTS ANALYSIS PIPELINE    
         ==============================================================================================
         Author                                                                      : ${params.author}
         Threshold AF                                                                : ${params.afmin}
         Threshold on known variants AF (gnomAD, EXAC, AVGP, SAHGP, TRYPANOGEN, KG)  : ${params.gnomad}
         Input dir                                                                   : ${in_dir}
         Results output                                                              : ${out_dir}
         ==============================================================================================
         
         """
.stripIndent()
println "============================================"
println "Preparing output folders for each population"
println "============================================"
POPS.each {
    out_it = file("${out_dir}/${it}")
	result = out_it.mkdirs()
	println result ? "Created ${it}" : "Cannot create folder: ${it}"
	out_it_stats = file("${out_dir_stats}/${it}")
	result = out_it_stats.mkdirs()
	println result ? "Created ${it}" : "Cannot create folder: ${it}"
	out_it_freq = file("${out_dir_freq}/${it}")
	result = out_it_freq.mkdirs()
	println result ? "Created ${it}" : "Cannot create folder: ${it}"
	out_it_impact = file("${out_dir_impact}/${it}")
	result = out_it_impact.mkdirs()
	println result ? "Created ${it}" : "Cannot create folder: ${it}"
}
println "=================================================================="
println "Listing vcfs files for each populations pattern ${params.alltype}"
println "=================================================================="
POPS.each {
    println "VCF Files in ${it} :"
	listWithHidden = file("${in_dir}/${it}/*_Eagle.baylor.*${params.alltype}.vcf.gz", hidden: false)
	listWithHidden.each {
		vcfs.add("$it")
   		println "File ${it}"
	}
	mappopvcfs.put("${it}", listWithHidden);
	println "Freq Count Files in ${it} :"
	listWithHidden = file("${in_dir_freq}/${it}/*_Eagle.baylor.*${params.alltype}.daf.frq.count", hidden: false)
	listWithHidden.each {
		ffreqcounts.add("$it")
   		println "File ${it}"
	}
	mappopffreqcounts.put("${it}", listWithHidden);
}
vcfs_ch = Channel.from(vcfs)
ffreqcounts_ch = Channel.from(ffreqcounts)
println "=================================================================="
println "Annotating Novel Variants in vcfs files for each populations         "
println "=================================================================="
process annNovelVariants {
    tag { "Annotating Novel Variants" } 
  
    input:
		val(vcf) from vcfs_ch 

    output:
    	
		file "${file(vcf).getName()}" into vcfAnnNovel
		file "${file(vcf).getName()}" into vcfAnnNoveltwo
	script:
	 
	"""
	python ${params.scripts}/python/NovelAnnoVCF.py  --vcf-input ${vcf} --af-min ${params.afmin} --gnom-ad ${params.gnomad}  --output-folder ${params.output_dir}/NovelVariants --out ${file(vcf).getName()}
	
    """
}
println "=================================================================="
println "Statistics for Novel Variants                                     "
println "=================================================================="
process statsNovelVariants {
    tag { "Stats on  Novel Variants" } 
    input:
		file anvcf from vcfAnnNovel 

    output:
    	
		file "stats.csv" into csvstats
	script:
	 
	"""
	python ${params.scripts}/python/NovelStats.py  --vcf-input ${anvcf}  --output-folder ${params.output_dir}/NovelStats --out stats.csv
	
    """
}
println "==================================================================="
println "Statistics on Rare Variants, Singletons, doubletons and Tripletons "
println "==================================================================="
process statsRareSDTVariants {
    tag { "Stats Rare Variants SDT" } 
    //storeDir "${params.homeagd}/sdt"
    input:
		val(ffreqcount) from ffreqcounts_ch 

    output:
    	
		file "${file(ffreqcount).getName()}" into csvStatsSTD
		file "${file(ffreqcount).getName()}" into csvStatsSTDImpact
		//file "${file(ffreqcount).toAbsolutePath()}" into csvStatsSTDPaths
	script:
	 
	"""
	python  ${params.scripts}/python/RareVariantsSDTtons.py  --freq-input ${ffreqcount}  --output-folder ${params.output_dir}/VariantsFreq --out ${file(ffreqcount).getName()}
	
    """
}
println "==========================================================================="
println "Combine Statistics on Rare Variants, Singletons, doubletons and Tripletons "
println "==========================================================================="
process combineStatsRareSDTVariants {
    tag { "Combining Stats Rare Variants SDT" } 
    
    input:
		file statsfreqs from csvStatsSTD.toList() 

    output:
		stdout dehors    	
	script:
		"""
		python  ${params.scripts}/python/CombineFreqs.py  --freqs-folder ${params.output_dir}/VariantsFreq 
		
		""" 
	
}
dehors.subscribe{println it}
println "==========================================================================="
println "Effect Impacts for Rare Variants, Singletons, doubletons and Tripletons    "
println "==========================================================================="
process variantsImpacts {
    tag { "Combining Stats Rare Variants SDT" } 
    
    input:
		file statsfreqs from csvStatsSTDImpact.toList() 
		file annotatednovel from vcfAnnNoveltwo.toList()

    output:
		stdout dehorsout    	
	script:
		"""
		python  ${params.scripts}/python/RareVariantsImpactBatch.py  --freqs-folder ${params.output_dir}/VariantsFreq --vcfs-folder ${params.output_dir}/NovelVariants --python-script ${params.scripts}/python/RareVariantsImpact.py --output-folder ${params.output_dir}/VariantsImpact
		
		""" 
	
}
dehorsout.subscribe{println it}
workflow.onComplete { 
	println ( workflow.success ? "Done annotating novel variants!" : "Something went south" )
}
