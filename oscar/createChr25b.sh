#Comment out the first two (Mali and Zambia below)as we already have them in createChr25.sh
#
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


#####Expanding to the rest
#Find all Benin freq files and concatenate them, removing the first line
find . -name "Benin*" | xargs -n 1 tail -n +2 > 2Benin_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
#Make column 1, chromosome into 25
awk '$1=25' 2Benin_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Benin_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Find all Botswana freq files and concatenate them, removing the first line
find . -name "Botswana*" | xargs -n 1 tail -n +2 > 2Botswana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq

#Make column 1, chromosome into 25
awk '$1=25' 2Botswana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Botswana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Deleting the leading 2 in the file to allow the analysis code to work
mv 2Benin_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Benin_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
mv 2Botswana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Botswana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq

#Find all Burkina freq files and concatenate them, removing the first line
find . -name "Burkina*" | xargs -n 1 tail -n +2 > 2Burkina_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
#Make column 1, chromosome into 25
awk '$1=25' 2Burkina_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Burkina_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Find all Cameroon freq files and concatenate them, removing the first line
find . -name "Cameroon*" | xargs -n 1 tail -n +2 > 2Cameroon_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq

#Make column 1, chromosome into 25
awk '$1=25' 2Cameroon_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Cameroon_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Deleting the leading 2 in the file to allow the analysis code to work
mv 2Burkina_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Burkina_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
mv 2Cameroon_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Cameroon_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq

#Find all Ghana freq files and concatenate them, removing the first line
find . -name "Ghana*" | xargs -n 1 tail -n +2 > 2Ghana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
#Make column 1, chromosome into 25
awk '$1=25' 2Ghana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Ghana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Find all Nigeria freq files and concatenate them, removing the first line
find . -name "Nigeria*" | xargs -n 1 tail -n +2 > 2Nigeria_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq

#Make column 1, chromosome into 25
awk '$1=25' 2Nigeria_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq > Nigeria_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.txt

#Deleting the leading 2 in the file to allow the analysis code to work
mv 2Nigeria_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Nigeria_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq
mv 2Ghana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq Ghana_phased_chr25_dp6_anc_f_dbsnp_snpeff.daf.frq