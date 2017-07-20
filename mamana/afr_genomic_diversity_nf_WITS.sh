#!/bin/bash

. $HOME/.bashrc

HOMEDIR="/home/mamana/GAPW/afr_genomic_diversity"
OUTDIR="/spaces/mamana/GAPW"

## load Python virtual environment
source activate ngs

## Nextflowscript here
cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/afr_genomic_diversity.nextflow.log \
    run ${HOMEDIR}/afr_genomic_diversity.nf \
    -c ${HOMEDIR}/afr_genomic_diversity_nf_WITS.config \
    -w ${OUTDIR}/work \
    -resume \
    -profile pbs
