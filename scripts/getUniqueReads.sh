#! /bin/bash
if [ "$#" -ne 3 ]; then
        echo "Illegal number of parameters"
        echo "filterToBam.sh <directory containing directories that contain sam/bam files> <output output directory> <email>"
        exit 1
fi

inputDir=$1
outputDir=$2
email=$3

# make output directory if it doesn't exist
if [ ! -d $outputDir ]
then
	mkdir $outputDir
fi


# loop for sam files
find $inputDir -name '*.sam' | while read samFile; do
	sampleName=${samFile/.sam/}
	sampleName=${sampleName##*/}

	echo "samtools view -Sbq 1 $samFile > $outputDir/${sampleName}.unique.bam"
	command="samtools view -Sbq 1 $samFile > $outputDir/${sampleName}.unique.bam"
	        # create qsub file
        echo "#!/bin/bash
#PBS -q hotel
#PBS -N $sampleName
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -o ${sampleName}_torque_output.txt
#PBS -e ${sampleName}_torque_error.txt
#PBS -V
#PBS -M $email
#PBS -m abe
#PBS -A glass-group
$command" > $outputDir/${sampleName}.torque.sh

qsub $outputDir/${sampleName}.torque.sh

done

# loop for bam files
find $inputDir -name '*.bam' | while read bamFile; do
	sampleName=${bamFile/.bam/}
	sampleName=${sampleName##*/}

	echo "samtools view -bq 1 $bamFile > $outputDir/${sampleName}.unique.bam"
	command="samtools view -bq 1 $bamFile > $outputDir/${sampleName}.unique.bam"
        # create qsub file
        echo "#!/bin/bash
#PBS -q hotel
#PBS -N $sampleName
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -o ${sampleName}_torque_output.txt
#PBS -e ${sampleName}_torque_error.txt
#PBS -V
#PBS -M $email
#PBS -m abe
#PBS -A glass-group
$command" > $outputDir/${sampleName}.torque.sh

qsub $outputDir/${sampleName}.torque.sh

done
