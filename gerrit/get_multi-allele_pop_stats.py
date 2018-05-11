#!/usr/bin/env python3
# This program is very specific to what we need to do
# 1) Generates allele counts
# ---- Allele counts ----
# Population	Bi-allelic count	Tri-allelic count	Other-allelic count	Allele types
# 2) Generate allele scenario counts
# Population	Two common alleles and one rare allele	Two rare alleles and one common allele

import vcf
import subprocess
import csv

def main():

#    base_path = '/home/gerrit/h3a/h3africa/gapw/african_diversity/manuscript_preperation/multi-allele/set1'
    base_path = '/spaces/gapw/diversity/gerrit/no_indels_combined_per_pop_vcfs'
    populations = {'BBC','BOT','BRN','BSZ','FNB','MAL','WGR'}
    population_info = {'BBC' : {},'BOT' : {},'BRN' : {},'BSZ' : {},'FNB' : {},'MAL' : {},'WGR' : {}}

    population_info = get_population_info_dictionary(base_path,populations,population_info)

    calculate_allele_counts_per_pop(population_info)
    get_allele_scenarios(population_info)

def get_population_info_dictionary(base_path, populations, population_info):

    for population in populations:
        #file_path = base_path + "/" + population + '.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.no_ref_0.000.exon.tsv'
        file_path = base_path + "/" + population + '.all.baylor_post_vqsr_clean.multi_alleles_only.no_indels.filterred_high_missing.related_removed.geDP6.no_ref_0.000.tsv'
        csv_file =  open(file_path,newline='')
        csv_reader = csv.reader(csv_file, delimiter='\t')
        csv_reader.__next__() # Skip header
        allele_info = {};
        for row in csv_reader:
            # Get dict key, the chromosome coord
            coord = ""
            coord = (row[0])
            ## Get alleles
            alleles = []
            alleles = row[1].split(',')
            ## Get allele frequences
            allele_freqs = []
            allele_freqs = row[2].split(',')
            # Append allle fequencies o alleles and add to allele info dict
            # Will look something like this for tri-alleles
            # 22:51207950 ['T', 'G', 'C', '0.982', '0.000', '0.018']
            allele_info[coord] = alleles + allele_freqs
        # Population population dict with allle info retrieved
        population_info[population] = allele_info

    return (population_info)

def calculate_allele_counts_per_pop(population_info):
    print ("---- Allele counts ----");
    print ("Population" + "\t" + "Bi-allelic count" + "\t" + "Tri-allelic count" + "\t" + "Other-allelic count" + "\t" + "Allele types")
    allele_types = []
    for population in list(population_info.keys()) :
        bi_allelic_count = 0
        tri_allelic_count = 0
        other_allelic_count = 0
        for coord,allele_info in population_info[population].items():
            other_allelic_count += 1
            # Just keep track of what allele typs we have e.g. bi,tri,quadri
            if not (round(len(allele_info) / 2) in allele_types):
                allele_types.append(round(len(allele_info) / 2))
            # Only focus on tri-alleles for now
            if((len(allele_info) / 2) == 3.0):
                tri_allelic_count += 1
                other_allelic_count -= 1
                # Freqs are in this range for tri-alleles
                for i in range (3,6):
                    if (float(allele_info[i]) == 0.000):
                        bi_allelic_count += 1
                        tri_allelic_count -= 1
                        if not (2 in allele_types):
                            allele_types.append(2)
        allele_types.sort()
        print (population + "\t" + str(bi_allelic_count) + "\t" + str(tri_allelic_count) + "\t" + str(other_allelic_count) + "\t%s" % ",".join(map(str, allele_types)))
    print ("--------------------------")

def get_allele_scenarios(population_info):
    print ("---- Allele scenarios counts ----");
    print ("Population" + "\t" + "Two common alleles and one rare allele" + "\t" + "Two rare alleles and one common allele")
    for population in list(population_info.keys()) :
        scenario_two_common_one_rare = 0
        scenario_one_common_two_rare = 0
        for coord,allele_info in population_info[population].items():
            # Only focus on tri-alleles for now
            if((len(allele_info) / 2) == 3.0):
                common_count = 0
                rare_count = 0
                for i in range (3,6):
                    if (float(allele_info[i]) >= 0.05):
                        common_count += 1
                    if (round(float(allele_info[i]),2) == 0.01):
                        rare_count += 1
                if (common_count == 2 and rare_count == 1):
                    # print("Two common and one rare" + "\t" + coord + "\t%s" % ",".join(map(str, allele_info)))
                    scenario_two_common_one_rare += 1
                if (common_count == 1 and rare_count == 2):
                    # print("One common and two rare" + "\t" + coord + "\t%s" % ",".join(map(str, allele_info)))
                    scenario_one_common_two_rare += 1

        print (population + "\t" + str(scenario_two_common_one_rare) + "\t" + str(scenario_one_common_two_rare))
    print ("-------------------------------------")

if __name__ == "__main__":
    main()
