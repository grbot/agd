#!/usr/bin/env python3
# Reads a VCF or GZipped VCF and reformats it into
# Chromosome position, alleles, allele frequences and sample names
# Output is printed to stdout
# E.g.
#Chrom:Pos  Alleles  Freq               430-ZLIU  431-ZLIU  432-ZLIU  433-ZLIU  434-ZLIU ...
#1:11232    C,A,T    0.958,0.027,0.015  0/0       0/0       0/0       0/0       0/0      ...
#2:15213    A,C,T    0.977,0.022,0.002  0/0       0/0       0/0       0/0       0/0      ...
#3:15555    A,T,C    0.052,0.404,0.544  2/2       1/2       1/2       1/2       2/2      ...

import vcf
from optparse import OptionParser

def main():
    parser = OptionParser(usage="usage: %prog -v input.vcf / input.vcf.gz")
    parser.add_option("-v", "--vcf", dest="vcf_file", help="Path to VCF to be reformatted.")
    (options, args) = parser.parse_args()
    if not options.vcf_file:
        print("No VCF specified, please specify with the -v flag.")
        return -1

    vcf_file = options.vcf_file
    vcf_reader = vcf.Reader(filename=vcf_file)

    header_printed = 0;

    for record in vcf_reader:

        if (header_printed == 0):
            # Print header chromosome position, alleles, allele frequences and sample names
            print ("Chrom:Pos\tAlleles\tFreq\t", end="\t")
            for i in range(0,len(record.samples) - 1):
                print (record.samples[i].sample, end='\t')
            print (record.samples[len(record.samples) - 1].sample)
            header_printed = 1

        # Print site data
        ## Print chromosome position
        print (record.CHROM, end=':')
        print (record.POS, end='\t')

        ## Print alleles (starting with reference)
        print (record.REF, end=',')
        for i in range(0,len(record.ALT) - 1):
            print (record.ALT[i], end=',')
        print (record.ALT[len(record.ALT) - 1], end='\t')

        ## Print allele frequencies but recalculate them first
        an = record.INFO['AN']
        ac_ref = (an - sum(record.INFO['AC'])) / an
        print ("{:.3f}".format(ac_ref), end=',')
        for i in range(0,len(record.INFO['AC']) - 1):
            print ("{:.3f}".format(record.INFO['AC'][i] / an), end=',')
        print ("{:.3f}".format(record.INFO['AC'][len(record.INFO['AC']) - 1] / an), end='\t')

        ## Print genotype data
        for i in range(0,len(record.samples) - 1):
            print (record.samples[i]['GT'], end='\t')
        print (record.samples[len(record.samples) - 1]['GT'])

if __name__ == "__main__":
    main()
