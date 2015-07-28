#! /bin/bash

# parse the input
inputDir=$1

if [ "$#" -ne 1 ]; then
	echo "Illegal number of parameters"
	echo "calcPBC.sh <directory containing samtool pileups> "
	exit 1
fi

find $inputDir -name '*.pileup' | while read pileupFile; do
	sampleName=${pileupFile%%.*}
	sampleName=${sampleName##*/}
	PBC=$(awk 'BEGIN {N1=0;ND=0} {if($4==1){N1+=1} ND+=1} END{print N1/ND}' ${pileupFile})
	echo -e "$sampleName\t$PBC"
done
