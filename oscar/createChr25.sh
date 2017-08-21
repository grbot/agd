#Find all Mali freq files and concatenate them, removing the first line
find . -name "Mali*" | xargs -n 1 tail -n +2 > 2Mali_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
#Make column 1, chromosome into 25
awk '$1=25' 2Mali_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Mali_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Find all Zambia freq files and concatenate them, removing the first line
find . -name "Zambia*" | xargs -n 1 tail -n +2 > 2Zambia_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq

#Make column 1, chromosome into 25
awk '$1=25' 2Zambia_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Zambia_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Deleting the leading 2 in the file to allow the analysis code to work
mv 2Zambia_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Zambia_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
mv 2Mali_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Mali_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq