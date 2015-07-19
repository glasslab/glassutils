#! /bin/bash

# parse the input
inputDir=$1

if [ "$#" -ne 1 ]; then
	echo "Illegal number of parameters"
	echo "calcPBC.sh <directory containing samtool pileups> "
	exit 1
fi

find $inputDir -name '*.sam' | while read samFile; do
	echo $samFile
	sampleName=${samFile%%.*sam}
	command=""
	if [ "$dataType" = "chip" ]
	then
		command="makeTagDirectory $sampleName -genome $genome -checkGC $samFile -format sam"

	elif [ "$dataType" = "rna" ]
	then
		command="makeTagDirectory $sampleName -genome $genome -checkGC $samFile -format sam -flip"
	fi
	# move tag directory to $input/tagDirs
	command="$command;mv $sampleName/ $inputDir/tagDirs/"

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
$command" > ${sampleName}.torque.sh

qsub ${sampleName}.torque.sh
done
