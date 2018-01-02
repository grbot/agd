#!/bin/bash

nextflow -log variants_analysis.nextflow.log \
    run variants_analysis.nf \
    -c variants_analysis.config \
    -w work \
    -resume \
-profile pbs
