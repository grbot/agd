echo "python $1 --vcf-input $2 --af-min 0.0 --gnom-ad 0.0 --output-folder $3"  | qsub -N $6 -l nodes=1:ppn=1,walltime=2:00:00,mem=5GB -e $4 -o $5
