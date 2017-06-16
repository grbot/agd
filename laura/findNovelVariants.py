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
        print '-p: <pathLoc> specify path (path/) where vcf file is located'
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
for chrom in range(1,23): 
    sys.stderr.write('starting chromosome %d\n' % chrom)
    vcfFName = '%sbaylor_phased_chr%d_dp6_anc_f_dbsnp_snpeff.vcf.gz' % (pathLoc, chrom)
    vcfReader = vcf.Reader(filename=vcfFName)
    for variant in vcfReader:
        c = c + 1
        if variant.INFO['snp138'][0] == None:
            if len(variant.INFO['AF']) > 0:
                sys.stderr.write('chr %d position %d VariantType=%s\n' % (chrom, variant.POS, variant.INFO['VariantType'])
                sys.exit()
            novel = novel + 1
            outFile.write('%s\t%d\t%f\n' % (variant.CHROM, variant.POS, variant.INFO['AF'][0]) )
            afSNPnovel.append(variant.INFO['AF'][0])
        else:
            afSNPrsID.append(variant.INFO['AF'][0])
outFile.close()
sys.stderr.write('Found %d novel variants out of %d (%f)\n' % (novel, c, float(novel)/d ))
sys.stderr.write('Average frequency novel variants: %f\n' % float(sum(afSNPnovel))/len(afSNPnovel) )
sys.stderr.write('Average frequency of discovered variants: %f\n' % float(sum(afSNPrsID))/len(afSNPrsID) ) 
