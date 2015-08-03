#! /bin/bash

# parse the input
inputDir=$1
outputDir=$2
email=$3

# calculates the PBC coefficient for a set of sam files contained in the same directory

if [ "$#" -ne 3 ]; then
	echo "Illegal number of parameters"
	echo "calcPBC.sh <directory that contains sam/bam files> <outputDirectory> <email>"
	exit 2
fi

# make directory to store  output
if [ ! -d $outputDir ]
then
	mkdir $outputDir
fi

find $inputDir -name '*.sam' | while read samFile; do
	sampleName=${samFile/.fastq.mm10.bowtie2.sam/}
	sampleName=${sampleName##*/}

	uniqueFile=$outputDir/${sampleName}.unique.bam
	sortedFile=$outputDir/${sampleName}.sorted
	pileupFile=$outputDir/${sampleName}.pileup

	command="samtools view -Sbq 1 $samFile > $outputDir/${sampleName}.unique.bam;"
	command+="samtools sort $uniqueFile $sortedFile;"
	command+="samtools mpileup ${sortedFile}.bam > $pileupFile;"
        command+="PBC=\$(awk 'BEGIN {N1=0;ND=0} {if(\$4==1){N1+=1} ND+=1} END{print N1/ND}' ${pileupFile});"
        command+="echo -e \"$sampleName\t\$PBC\" >>$outputDir/pbc.tsv;"

	echo "submitting $command"

	# create qsub file
	echo "#!/bin/bash
#PBS -q hotel
#PBS -N ${sampleName}_pbc
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -o $outputDir/${sampleName}_pbc_torque_output.txt
#PBS -e $outputDir/${sampleName}_pbc_torque_error.txt
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

	uniqueFile=$outputDir/${sampleName}.unique.bam
	sortedFile=$outputDir/${sampleName}.sorted
	pileupFile=$outputDir/${sampleName}.pileup

	command="samtools view -bq 1 $bamFile > $outputDir/${sampleName}.unique.bam;"
	command+="samtools sort $uniqueFile $sortedFile;"
	command+="samtools mpileup ${sortedFile}.bam > $pileupFile;"
        command+="PBC=\$(awk 'BEGIN {N1=0;ND=0} {if(\$4==1){N1+=1} ND+=1} END{print N1/ND}' ${pileupFile});"
        command+="echo -e \"$sampleName\t\$PBC\" >>$outputDir/pbc.tsv;"
	command+="if [ \$(ls -l $outputDir/*.pileup) -eq \$(wc -l $outputDir/pbc.tsv) ] then echo \"moving PBC results\"; fi"
	echo "submitting $command"

	# create qsub file
	echo "#!/bin/bash
#PBS -q hotel
#PBS -N ${sampleName}_pbc
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -o $outputDir/${sampleName}_pbc_torque_output.txt
#PBS -e $outputDir/${sampleName}_pbc_torque_error.txt
#PBS -V
#PBS -M $email
#PBS -m abe
#PBS -A glass-group
$command" > $outputDir/${sampleName}.torque.sh

qsub $outputDir/${sampleName}.torque.sh
done
