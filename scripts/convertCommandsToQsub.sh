#!/bin/bash

#converts individual lines of a script into individual qsub scripts

# parse the input
scriptFile=$1
outputDir=$2
wallTime=$3
email=$4

if [ "$#" -ne 4 ]; then
        echo "Illegal number of parameters"
        echo "create_tagDir.sh <original script file> <outputDirectory> <wallTime> <email>"
        exit 1
fi

# create output directory if it doesn't exist
if [ ! -d $outputDir ]
then
	mkdir $outputDir
fi
counter=1
while read line 
do
        echo "submitting $line"
        # create qsub file
        echo "#!/bin/bash
#PBS -q hotel
#PBS -N $counter
#PBS -l nodes=1:ppn=8
#PBS -l walltime=${wallTime}:00:00
#PBS -o ${counter}_torque_output.txt
#PBS -e ${counter}_torque_error.txt
#PBS -V
#PBS -M $email
#PBS -m abe
#PBS -A glass-group

	$line" > $outputDir/${counter}.torque.sh

	qsub $outputDir/${counter}.torque.sh
	counter=$((counter+1))

done < $scriptFile

