#!/usr/bin/perl -w

#my $bowtieDir = "/bioinformatics/software/bowtie/bowtie-0.12.7/";
my $bowtieDir = "/projects/glass-lab/bioinformatics/bowtie/bowtie-0.12.7/";
my $gtfDir = "/projects/glass-lab/tophat/gtfs/";
#my $gtfDir = "/Volumes/Unknowme/mapping/gtf/";

#my $libType = "fr-secondstrand"; # works best for our type of RNA-seq
my $libType = "fr-firststrand";

my $nice = 10;

# color space
#/bioinformatics/software/bowtie/bowtie-0.12.7/bowtie -e 200 --sam -m 1 -n 3 --chunkmbs 128 -q --best -p 6 -C /bioinformatics/software/bowtie/bowtie-0.12.7/indexes/mm9_c SRR039863.fastq  > SRR039863.fastq.sam


my $qualOptions = "";
#$qualOptions = " -e 200 ";
$qualOptions = " -e 200 --nomaqround";
#$qualOptions = " -e 100 --solexa1.3-quals ";

my $maxTopHatHits = 5;
my $fastaFlag = 0;
my $reduce = 0;
my $minLength = 14;
my $unmapFlag = 0;


if (@ARGV < 3) { 
	print STDERR "\n\tUsage: map-bowtie.pl <# cpu> <bowtie genome> [options] <FASTQ file 1> [FASTQ file 2] ...\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-max <#> (Maximum read length, no max)\n";
	print STDERR "\t\t-min <#> (Minimum read length, default: $minLength)\n";
	print STDERR "\t\t-f (input sequences are FASTA files)\n";
	print STDERR "\t\t-un (output unaligned sequences)\n";
	print STDERR "\t\t-tophat (Perform spliced alignment with TopHat)\n";
	print STDERR "\t\t\t--segment-length <#> (set to half the read length if less than 50, default 25)\n";
	print STDERR "\t\t-junc <file> (Input splice junction file for TopHat [also turns off novel-juncs])\n";
	print STDERR "\t\t-reduce <#> (Keep reducing sequence length by # until sequence maps, assumes 50bp max)\n";
	print STDERR "\t\t-uniq (report only unique alignment, default)\n";
	print STDERR "\t\t-top2 (report top 2 sequences in bowtie format)\n";
	print STDERR "\t\t-sam (report alignment in SAM format, default)\n";
	print STDERR "\t\t-bowtie (report bowtie format)\n";
	print STDERR "\t\t-C (color space)\n";
	print STDERR "\t\t-nice <#> (nice priority level, -20 highest priority, 20 lowest, default: 10\n";
	print STDERR "\n\tAvaliable Genomes:\n";

	`ls "$bowtieDir"/indexes/*rev.1.ebwt > .ls`;
	open IN, ".ls";
	while (<IN>) {
		chomp;
		$_ =~ s/^.*indexes\///;
		$_ =~ s/\.rev\.1\.ebwt//;
		print "\t\t$_\n";
	}
	close IN;
	`rm .ls`;	
	print "\n";
	exit;
}

my $maxMismatches = 3;
my $cpu = $ARGV[0];
my $genome = $ARGV[1];
my $len = 0;
my $colorSpace = '';
my $tophatFlag = 0;
my $juncFile = "";
my $format = "--sam";
my $top2Flag = "-m 1";
my $tophatSegmentLength = 25;
my $tophatSegmentMisMatch = 2;


my @seqFiles = ();
for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-max') {
		$len = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minLength = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-f') {
		$fastaFlag = 1;
	} elsif ($ARGV[$i] eq '-C') {
		$colorSpace = "-C";
	} elsif ($ARGV[$i] eq '-junc') {
		$juncFile .= " -j $ARGV[++$i] --no-novel-juncs ";
		next;
	} elsif ($ARGV[$i] eq '-nice') {
		$nice = $ARGV[++$i];
		print STDERR "\tNice level: $nice\n";
	} elsif ($ARGV[$i] eq '-un') {
		$unmapFlag = 1;
	} elsif ($ARGV[$i] eq '-reduce') {
		$reduce = $ARGV[++$i];
		if ($len < 1) {
			$len = 50;
		}
	} elsif ($ARGV[$i] eq '-uniq') {
		$top2Flag = "-m 1";
	} elsif ($ARGV[$i] eq '--segment-length') {
		$tophatSegmentLength = $ARGV[++$i];
		if ($tophatSegmentLength < 25) {
			$tophatSegmentMisMatch = 1;
		}
	} elsif ($ARGV[$i] eq '-top2') {
		$top2Flag = "-k 2";
		$format = "";
	} elsif ($ARGV[$i] eq '-bowtie') {
		$format = "";
	} elsif ($ARGV[$i] eq '-sam') {
		$format = "--sam";
	} elsif ($ARGV[$i] eq '-tophat') {
		$tophatFlag = 1;
		if ($minLength < 20) {
			$minLength = 20;
		}
		next;
	} else {
		print STDERR "\tWill align $ARGV[$i]\n";
		push(@seqFiles, $ARGV[$i]);
	}
}

print STDERR "\n\t-----------------------------------------\n";
print STDERR "\tAligning to $genome with $cpu CPUs\n\n";


$program = $bowtieDir  .  "/bowtie";
#print STDERR "Current Program: $program\n";
#print STDERR "Will map $len bp\n";

my $tophatTmpDir = rand() . ".tophat";

my @pids = ();
foreach (@seqFiles) {

	my $inputfile = $_;
	my $zipFlag = 0;
	if ($inputfile =~ /\.gz$/) {
		print STDERR "\tUnzipping file $inputfile\n";
		`gunzip "$inputfile"`;
		$zipFlag = 1;
		$inputfile =~ s/\.gz$//;
	}

	my $outputfile = $inputfile. ".$genome";
	$outputfile .= ".$len" if ($len > 0);
	if ($format eq '--sam') {
		$outputfile .= ".sam";
	} else {
		$outputfile .= ".align";
	}

	my $lenFile = $inputfile . ".all.bp";
	my $lenbp = "all.bp";
	if ($len > 0) {
		$lenFile = $inputfile . ".$len.bp";
		$lenbp = $len . ".bp";
	}
	my $lenFile2 = $inputfile . ".lengths";

	my $unmapOptions = '';
	if ($unmapFlag) {
		if ($len > 0) {
			$unmapOptions = " --un " . $inputfile . ".$genome.$len.unalign.fq";
		} else {
			$unmapOptions = " --un " . $inputfile . ".$genome.unalign.fq";
		}
	}

	# trim sequences to get them the right size
	my $opt = "";
	$opt = "-len $len" if ($len > 0);
	`nice -n $nice homerTools trim $opt -min $minLength -suffix $lenbp $inputfile`;

	if ($tophatFlag == 1) {
		$outputfile = $inputfile . ".$genome.tophat.bam";
		$outputfile2 = $inputfile . ".$genome.tophat.junc";
		`nice -n $nice tophat -o "$tophatTmpDir" -g $maxTopHatHits -G $gtfDir/$genome.refseq.gtf $juncFile -F 0 --library-type $libType --segment-mismatches $tophatSegmentMisMatch --segment-length $tophatSegmentLength -p $cpu $bowtieDir/indexes/$genome $lenFile`;
		`mv "$tophatTmpDir/accepted_hits.bam" "$outputfile"`;
		`parseTophatJunctions.pl "$tophatTmpDir/junctions.bed" > "$outputfile2"`;
		`rm -r "$tophatTmpDir"`;
	} else {
		if ($reduce == 0) {
			if ($fastaFlag == 1) {
				`nice -n $nice $program $colorSpace -n $maxMismatches $unmapOptions $format --chunkmbs 128 -f --best $top2Flag -p $cpu $genome $lenFile $outputfile 2> $outputfile.mapstats`;	
			} else {
				#for fastq
				`nice -n $nice $program $colorSpace -n $maxMismatches $unmapOptions $format --chunkmbs 128 -q --best $top2Flag $qualOptions -p $cpu $genome $lenFile $outputfile 2> $outputfile.mapstats`;	
			}
		} else {


			#print STDERR "`mv $lenFile tmp.unalign`;\n";
			my $tmpUnalign = "$outputfile.tmp.unalign";
			my $tmpIn = "$outputfile.tmp.in";
			
			`mv $lenFile $tmpUnalign`;
			my $maxReduce = $len - $minLength;
			for (my $j=0;$j<=$maxReduce;$j+=$reduce) {

				`mv $tmpUnalign $tmpIn`;

				#`homerTools trim -len $j -min $minLength -suffix $j $tmpUnalign`;

				$unmapOptions = " --un $tmpUnalign";

				if ($fastaFlag == 1) {
					`nice -n $nice $program $colorSpace -3 $j -n $maxMismatches $unmapOptions $format --chunkmbs 128 -f --best $top2Flag -p $cpu $genome $tmpIn $outputfile.trim$j 2> $outputfile.mapstats`;	
				} else {
					#for fastq
					`nice -n $nice $program $colorSpace -3 $j -n $maxMismatches $unmapOptions $format --chunkmbs 128 -q --best $top2Flag $qualOptions -p $cpu $genome $tmpIn $outputfile.trim$j 2> $outputfile.mapstats`;	
				}
			}
			`cat $outputfile.trim* > $outputfile`;
			my $XX = $inputfile . ".$genome.unalign.fq";
			`mv $tmpUnalign $XX`;


		}
		open IN, "$outputfile.mapstats";
		while (<IN>) {
			print STDERR "$_";
		}
		close IN;
		
	}
	`rm $lenFile $lenFile2`;
	if ($zipFlag) {
		`gzip $inputfile`;
	}
}


