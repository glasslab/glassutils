#! /bin/bash

# parse the input
inputDir=$1
dataType=$2
genome=$3
email=$4

if [ "$#" -ne 4 ]; then
	echo "Illegal number of parameters"
	echo "create_tagDir.sh <directory containing directories that contain sam files> <chip|rna> <genome> <email>"
	exit 1
fi

# make directory to store tag directories
if [ ! -d $inputDir/tagDirs ]
then
	mkdir $inputDir/tagDirs
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
