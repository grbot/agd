import getopt
import sys

import vcf

pathLoc = None
outFName = None

options, args = getopt.getopt(sys.argv[1:], 'p:o:h')

for opt, par in options:
    if opt == '-h':	
        print 'options:'	
        print '-h: print this help and exit'
        print '-p: <pathLoc> specify vcf file'
        print '-o: <outFName> specify output txt file'
        sys.exit()
    elif opt == '-p':
        pathLoc = par
    elif opt == '-o':
        outFName = par
    else:
        raise StandardError, 'unhandled option "%s"' % opt 

outFile = open(outFName, 'w')
outFile.write('CHR\tPOS\tFREQ\n')
c = 0 #counter of total number of polymorphisms
novel = 0 #counter of novel variants
afSNPrsID = []
afSNPnovel = []
for chrom in [21]: #range(1,23): 
    sys.stderr.write('starting chromosome %d\n' % chrom)
    vcfFName = '%sGhana_phased_chr%d_dp6_anc_f_dbsnp_snpeff.vcf.gz' % (pathLoc, chrom)
    vcfReader = vcf.Reader(filename=vcfFName)
    for variant in vcfReader:
        c = c + 1
        if variant.INFO['snp138'][0] == None: # Variant not in dbSNP and therefore, novel
            if len(variant.INFO['AF']) > 1:
                print variant.INFO
                print variant.ALT
                print variant.REF
                sys.stderr.write('chr %d position %d VariantType=%s\n' % (chrom, variant.POS, variant.INFO['VariantType']))
                sys.exit()
            novel = novel + 1
            outFile.write('%s\t%d\t%d\t%f\n' % (variant.CHROM, variant.POS-1, variant.POS, variant.INFO['AF'][0]) )
            afSNPnovel.append(float(variant.INFO['AF'][0]))
        else:
            afSNPrsID.append(float(variant.INFO['AF'][0]))
outFile.close()
print afSNPnovel[0:10]
sys.stderr.write('Found %d novel variants out of %d (percent: %f)\n' % (novel, c, float(novel)/c ))
sys.stderr.write('Average frequency novel variants: %f\n' % (sum(afSNPnovel)/len(afSNPnovel)) )
sys.stderr.write('Average frequency of known variants: %f\n' % (sum(afSNPrsID)/len(afSNPrsID)) ) 