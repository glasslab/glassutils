#!/bin/bash
inputDir=$1
for samFile in $inputDir/*sam
do
	bamFile=${samFile/.sam/.bam}
	echo "samtools view -bS $samFile > $bamFile"
	samtools view -bS $samFile > $bamFile
done

