#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=3
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -M mbymam001@myuct.ac.za
#PBS -m ea
#PBS -o /spaces/mamana/AFRICA_CHIP/GAPW/LOG/afr_genomic_diversity.out

. $HOME/.bashrc

HOMEDIR="${HOME}/GAPW/afr_genomic_diversity"
OUTDIR="/mnt/lustre/users/mmbiyavanga/GAPW"

## load Python virtual environment
source activate ngs_py35

## Nextflowscript here
cd ${HOMEDIR}
nextflow -log ${OUTDIR}/LOG/afr_genomic_diversity.nf.log \
    run ${HOMEDIR}/afr_genomic_diversity.nf \
    -c ${HOMEDIR}/afr_genomic_diversity.nf_CHPC.config \
    -w ${OUTDIR}/work \
    -resume
#    -profile pbs
