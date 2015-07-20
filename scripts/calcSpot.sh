#! /bin/bash

# parse the input
inputDir=$1
outputDir=$2
email=$3

if [ "$#" -ne 3 ]; then
	echo "Illegal number of parameters"
	echo "calculatePileups.sh <directory that contains sam/bam files> <directory that contains bed files> <outputDirectory> <email>"
	exit 1
fi

# make directory to store  output
if [ ! -d $outputDir ]
then
	mkdir $outputDir
fi

find $inputDir -name '*.sam' | while read samFile; do
	sampleName=${samFile/.fastq.mm10.bowtie2.sam/}
	sampleName=${sampleName##*/}

	sortedFile=$outputDir/${sampleName}.sorted

	pileupFile=$outputDir/${sampleName}.pileup

	command="samtools sort $samFile $sortedFile"
	command+="; samtools mpileup ${sortedFile}.bam > $pileupFile"

	echo "submitting $command"

	# create qsub file
	echo "#!/bin/bash
#PBS -q hotel
#PBS -N $sampleName
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -o ${sampleName}_torque_output.txt
#PBS -e ${sampleName}_torque_error.txt
#PBS -V
#PBS -M $email
#PBS -m abe
#PBS -A glass-group
$command" > $outputDir/${sampleName}.torque.sh

qsub $outputDir/${sampleName}.torque.sh
done


find $inputDir -name '*.bam' | while read bamFile; do
	sampleName=${bamFile/.fastq.mm10.bowtie2.bam/}
	sampleName=${sampleName##*/}

	sortedFile=$outputDir/${sampleName}.sorted

	pileupFile=$outputDir/${sampleName}.pileup

	command="samtools sort $bamFile $sortedFile"
	command+="; samtools mpileup ${sortedFile}.bam > $pileupFile"

	echo "submitting $command"

	# create qsub file
	echo "#!/bin/bash
#PBS -q hotel
#PBS -N $sampleName
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -o ${sampleName}_torque_output.txt
#PBS -e ${sampleName}_torque_error.txt
#PBS -V
#PBS -M $email
#PBS -m abe
#PBS -A glass-group
$command" > $outputDir/${sampleName}.torque.sh

qsub $outputDir/${sampleName}.torque.sh
done
