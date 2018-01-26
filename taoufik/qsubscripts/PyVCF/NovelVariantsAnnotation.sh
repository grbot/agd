#!/bin/bash

MainJob="Extracting Novel Variants By annotation"
echo $MainJob
annotatedvcfpath=$1
resultsoutput=$2
runhome=$(pwd)
echo "VCF path :"$annotatedvcfpath
echo "Results path :"$resultsoutput
if [ -d $resultsoutput ] ; then
	echo $resultsoutput" exists"
else
	mkdir $resultsoutput
fi
cd $annotatedvcfpath
for dir in *;
	do
	echo $dir
	cd $dir
	for fvcf in *;
		do
		echo "File VCF -- "$fvcf
		sh $runhome/NovelAnnoVCF.sh $runhome/NovelAnnoVCF.py $annotatedvcfpath$dir"/"$fvcf $resultsoutput $fvcf"-error.log" $fvcf"-out.log" "Job-Novel-Variants-in-"$fvcf
		# ""| qsub -N $fvcf -l nodes=1:ppn=1,walltime=2:00:00,mem=5GB -e $fvcf"-error.log" -o $fvcf"-out.log" "
	done
	cd ..
done
