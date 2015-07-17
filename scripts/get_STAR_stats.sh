inputPath=$1 # path to a directory containing direcotries that contain a STAR log file

echo -e "Sample Name\tNumber of input reads\tUniquely mapped reads %\t% of reads mapped to multiple loci"
for dir in $inputPath/*/
do
	sampleName=${dir%%\/}
	sampleName=${sampleName#*/}
	if [ -f $dir/Log.final.out ]
	then
		ir=$(awk '/Number of input reads/ {print $6}' $dir/Log.final.out) # number of input reads
		umr=$(awk '/Uniquely mapped reads \%/ {print $6}' $dir/Log.final.out) # uniquely mapped reads fraction
		mmr=$(awk '/\% of reads mapped to multiple loci/{print $9}' $dir/Log.final.out) # multimapped fraction
		toPrint="$sampleName\t$ir\t$umr\t$mmr"
		toPrint=${toPrint//\%/}
		echo -e $toPrint
	else
		# print to stderr
		>&2 echo "no log file found for $sampleName"
	fi
done
