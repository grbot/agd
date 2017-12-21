#!/bin/bash

nextflow -log  pipeline.log run variants_analysis.nf  -c variants_analysis.config --datavcf /spaces/gapw/diversity/mamana/VCF_POP/BAYLOR/VCF --datafreq /spaces/mamana/GAPW/VCF_POP/BAYLOR/DAF  -profile pbs
