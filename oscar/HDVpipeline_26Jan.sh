#!/bin/bash
#cd into analysis directory
cd /spaces/gapw/diversity/oscar/HDV

#runcolumn7.sh makes a 7the column and keeps there chromoosme information in columnn 1
#renamecol1.sh codes column 1 as chr25 to allow HDV_Figure1_VCFTools6test5.R  to give results for all chromosomes 1-22
#restorecol1.sh restores correct chromosome information back to column1

$pathtofrq_files='/spaces/gapw/diversity/mamana/VCF_POP/BAYLOR/DAF'
$pathtovcf_files='/spaces/gapw/diversity/mamana/VCF_FILTERED/BAYLOR'

declare -a pop =("BOT" "BRN" "CAM" "FNB" "WGR" "MAL" "ZAM")
  do
         for p in `seq 0 6`;
             do
             #makes a 7the column and keeps there chromoosme information in columnn 1
             awk '$7=$1' $pathtofrq_files/${pop[$p]}/${pop[$p]}_dp6_anc_f_dbsnp_snpeff.daf.frq > ${pop[$p]}_3dp6_anc_f_dbsnp_snpeff.daf.frq
             
             #renamecol1.sh codes column 1 as chr25 to allow HDV_Figure1_VCFTools6test5.R  to give results for all chromosomes 1-22
             awk '$1=25' ${pop[$p]}_3dp6_anc_f_dbsnp_snpeff.daf.frq > ${pop[$p]}_4dp6_anc_f_dbsnp_snpeff.daf.frq
             
             #restorecol1.sh restores correct chromosome information back to column1
             awk '$1=$7' ${pop[$p]}_4dp6_anc_f_dbsnp_snpeff.daf.frq > ${pop[$p]}_5dp6_anc_f_dbsnp_snpeff.daf.frq
             done
   done


#HDV_Figure1_VCFTools6test5.R Rscript which does the HDV calculations; generates 25HDVList_poppair.txt 
#and draws figures of HDV in PDF
Rscript HDV_Figure1_VCFTools6test5.R 25

#runhdv.sh correct chromosome results are in the fifth column; copy these into the first to replace chr25
# runhdv2.sh makes a tab delimited version of HDVlist; good for bedtools intersect
#HDV_bedtoolsintersect.sh  Intersect between the hdv file and 5kb windows keeping all the windows with HDVs
#sortbed.sh  sort output (seems to work and keep only highest per 5kb window)
#sorted_loci.sh Keeps only columns with HDV data: chr pos refallele altallele
#array_pos.sh; takes highest HDV within 5kb widnow ${arr[$k]}_HDV_5K2.bed from above and outputs 2 columns chr and position in ${arr[$k]}_HDV_5K2_pos.txt
declare -a poppair=("BOT-BRN" "BOT-CAM" "BOT-FNB" "BOT-MAL" "BOT-WGR" "BOT-ZAM" "BRN-CAM" "BRN-FNB" "BRN-MAL" "BRN-WGR" "BRN-ZAM" "CAM-FNB" "CAM-MAL" "CAM-WGR" "CAM-ZAM" "FNB-MAL" "FNB-WGR" "FNB-ZAM" "MAL-WGR" "MAL-ZAM" "WGR-ZAM")
declare -a arr=("BOT_BRN" "BOT_CAM" "BOT_FNB" "BOT_MAL" "BOT_WGR" "BOT_ZAM" "BRN_CAM" "BRN_FNB" "BRN_MAL" "BRN_WGR" "BRN_ZAM" "CAM_FNB" "CAM_MAL" "CAM_WGR" "CAM_ZAM" "FNB_MAL" "FNB_WGR" "FNB_ZAM" "MAL_WGR" "MAL_ZAM" "WGR_ZAM")

            for k in `seq 0 20`;
                do
                	#runhdv.sh correct chromosome results are in the fifth column; copy these into the first to replace chr25
                    awk '$1=$5' 25HDVList_${poppair[$k]}Delta60.txt > 25HDVList_${poppair[$k]}Delta60_2.txt   
                    
                    #runhdv2.sh make a tab delimited version of HDVlist; good for bedtools intersect
                    awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}'  25HDVList_${poppair[$k]}Delta60_2.txt  > 25HDVList_${poppair[$k]}Delta60_3.txt
                    
                    #HDV_bedtoolsintersect.sh  Intersect between the hdv file and 5kb windows keeping all the windows with HDVs
                    bedtools intersect -wa  -wb -a chrall_5K.bed  -b 	25HDVList_${poppair[$k]}Delta60_3.txt	>	${arr[$k]}_HDV_5K.bed
                    
                    #sortbed.sh  sort output check (seems to work and keep only highest per 5kb window)
                    sort -rnk3 ${arr[$k]}_HDV_5K.bed | awk 	'!x[$2]++'	>	${arr[$k]}_HDV_5K2.bed
                    #sorted_loci.sh Keeps only columns with HDV data: chr pos refallele altallele
                    awk '{print $4"\t"$5"\t"$7"\t"$8}'  ${arr[$k]}_HDV_5K2.bed  > ${arr[$k]}_HDV_sorted.txt
                    #array_pos.sh; takes highest HDV within 5kb widnow ${arr[$k]}_HDV_5K2.bed from above and outputs 2 columns chr and position in ${arr[$k]}_HDV_5K2_pos.txt
                    awk '{print $1"\t"$5}'  ${arr[$k]}_HDV_5K2.bed  >  ${arr[$k]}_HDV_5K2_pos.txt
             done


#Counting the number of lines per population pair before and after filtering with the 5kb filtering

#before5kb_filt.sh  count number of lines for all poplation pairs before 5kb filtering
#after5kb_filt.sh count number of lines for all population pairs after 5kb filtering
declare -a arr=("BOT_BRN" "BOT_CAM" "BOT_FNB" "BOT_MAL" "BOT_WGR" "BOT_ZAM" "BRN_CAM" "BRN_FNB" "BRN_MAL" "BRN_WGR" "BRN_ZAM" "CAM_FNB" "CAM_MAL" "CAM_WGR" "CAM_ZAM" "FNB_MAL" "FNB_WGR" "FNB_ZAM" "MAL_WGR" "MAL_ZAM" "WGR_ZAM")
for k in `seq 0 20`;
    do
      #before5kb_filt.sh  count number of lines for all poplation pairs before 5kb filtering
        wc -l ${arr[$k]}_HDV_5K.bed > before5kb_filt.txt

      #after5kb_filt.sh count number of lines for all population pairs after 5kb filtering
        wc -l ${arr[$k]}_HDV_5K2.bed > after5kb_filt.txt
    done

#Workflow_HDV_allpops.sh : Wrapper script for runvcfs.sh
#runvcfs.sh gets the filtered vcf positions file per population pair ${arr[$k]}_HDV_5K2_pos.txt and queries the vcf file to get just the infor in these positions

for i in `seq 1 22`; #i 
        #for i in `seq 15 15`; #i 
        do

                for k in `seq 0 20`;
                do
            
    vcftools \
    --gzvcf $pathtovcf_files/Eagle.baylor.${i}_snpeff_dbsnp_anc_gwascat_clinvar_cosmic_mafs-ExAC-KG-agvp-gnomAD-sahgp-trypan_noRelated_noHighMiss_dp6_ancOnly.vcf.gz \
    --positions ${arr[$k]}_HDV_5K2_pos.txt \
    --recode-INFO-all --recode --stdout \
    | gzip  > filter${arr[$k]}chr${i}.vcf.gz
                done

done

#The next steps required vcfquery ; which is part of the vcf.pm (perl module)
#This module failed to install on the UCT server
#I therefore 
#1. installed it on my laptop and ran it as follows:

#2.Copied the filtered vcf files into this directory
cd /Users/oscars_mac/Documents/GAPW/filtervcfs
scp oscar@cream.core.wits.ac.za:/home/oscar/HDV/filtervcfs.tar.gz .

#3.Extracted these vcfs within this directory
tar -xzvf filtervcfs2.tar.gz 

#4.Ran runquery_wrapper.sh; which calls the vcf perl module command

#command: sh /Users/oscars_mac/Documents/GAPW/Scripts/runquery_wrapper.sh

 #!/bin/bash
cd /Users/oscars_mac/Documents/GAPW/filtervcfs/
 declare -a arr=("BOT_BRN" "BOT_CAM" "BOT_FNB" "BOT_MAL" "BOT_WGR" "BOT_ZAM" "BRN_CAM" "BRN_FNB" "BRN_MAL" "BRN_WGR" "BRN_ZAM" "CAM_FNB" "CAM_MAL" "CAM_WGR" "CAM_ZAM" "FNB_MAL" "FNB_WGR" "FNB_ZAM" "MAL_WGR" "MAL_ZAM" "WGR_ZAM")

for i in `seq 1 22`; #i 
        #for i in `seq 15 15`; #i 
        do


                for k in `seq 0 20`;
                do
             #/Users/oscars_mac/Documents/GAPW/Scripts/runquery.sh $i $k
             gunzip filter${arr[$k]}chr${i}.vcf.gz
             /usr/local/bin/vcf-query -f '%CHROM %POS %ID  %INFO/ANN\n' filter${arr[$k]}chr${i}.vcf > ${arr[$k]}_output_${i}.txt

             done

done


#get_ensembl_ID.sh calls outfile3.sh above and makes a unique list of all genes per chromosome; 
#!/bin/bash
i=$1
k=$2
cd /Users/oscars_mac/Documents/GAPW/filtervcfs/
 declare -a arr=("BOT_BRN" "BOT_CAM" "BOT_FNB" "BOT_MAL" "BOT_WGR" "BOT_ZAM" "BRN_CAM" "BRN_FNB" "BRN_MAL" "BRN_WGR" "BRN_ZAM" "CAM_FNB" "CAM_MAL" "CAM_WGR" "CAM_ZAM" "FNB_MAL" "FNB_WGR" "FNB_ZAM" "MAL_WGR" "MAL_ZAM" "WGR_ZAM")

for i in `seq 1 22`; #i 
        #for i in `seq 15 15`; #i 
        do

                for k in `seq 0 20`;
                do
             sh /Users/oscars_mac/Documents/GAPW/Scripts/outfile3.sh $i $k > genelist_temp.txt
             sort genelist_temp.txt | uniq genelist_temp.txt>list_${arr[$k]}_chr${i}.txt
             done

done


#concat_HDVgenes.sh get all chromosomes with HDV genes identified into one file per population pair
 #!/bin/bash
 cd /Users/oscars_mac/Documents/GAPW/filtervcfs
 declare -a arr=("BOT_BRN" "BOT_CAM" "BOT_FNB" "BOT_MAL" "BOT_WGR" "BOT_ZAM" "BRN_CAM" "BRN_FNB" "BRN_MAL" "BRN_WGR" "BRN_ZAM" "CAM_FNB" "CAM_MAL" "CAM_WGR" "CAM_ZAM" "FNB_MAL" "FNB_WGR" "FNB_ZAM" "MAL_WGR" "MAL_ZAM" "WGR_ZAM")

            for k in `seq 0 20`;
                do
                        find . -name "list_${arr[$k]}*" | xargs -n 1 tail -n +2 > genesHDV_${arr[$k]}.txt
             done





#outfile3.sh called by get_ensembl_ID takes ${arr[$k]}_output_${i}.txt output from runquery.sh; reads it and outputs it to stdout for parsing to other commands
 #!/bin/bash
i=$1
k=$2
declare -a arr=("BOT_BRN" "BOT_CAM" "BOT_FNB" "BOT_MAL" "BOT_WGR" "BOT_ZAM" "BRN_CAM" "BRN_FNB" "BRN_MAL" "BRN_WGR" "BRN_ZAM" "CAM_FNB" "CAM_MAL" "CAM_WGR" "CAM_ZAM" "FNB_MAL" "FNB_WGR" "FNB_ZAM" "MAL_WGR" "MAL_ZAM" "WGR_ZAM")
cd /Users/oscars_mac/Documents/GAPW/filtervcfs/
file="/Users/oscars_mac/Documents/GAPW/filtervcfs/${arr[$k]}_output_${i}.txt"
while IFS= read -r line
do
        # display $line or do somthing with $line
#	printf '%s\n' "$line"
#done <"$file"


IFS='|' read -ra ANN <<< "$line"
#for i in "${ADDR[@]}"; do
    # process "$i"
#done

while IFS='|' read -ra ANN; do
      declare -a ANN="${ANN[@]}"
          #echo  "$i"
          echo ${ANN[4]}
      
 done <<< "$line"
done <"$file"