
my $method = 0;
#!/usr/bin/perl -w

use strict;
use File::Basename;

#map_on_tscc.sh STAR /data/archive/15-04-06-glass/Project_Eniko /oasis/tscc/scratch/vlink/eniko mm10 vlink vlink@ucsd.edu
#call this script
#perl create_tag_directories_on_tscc.pl STAR /data/archive/15-04-06-glass/Project_Eniko /oasis/tscc/scratch/vlink/eniko mm10 vlink vlink@ucsd.edu

if(@ARGV < 4) {
	print "Use the map_on_tscc_command options for this programm!\n";
	exit;
}
my $method = $ARGV[0];
my $path_glassome = $ARGV[1];
$path_glassome = substr($path_glassome, 1);
$path_glassome = "/projects/ps-glasslab-" . $path_glassome;
my $path_tscc = $ARGV[2];
my $genome = $ARGV[3];
my $rna_seq_old = 0;

if(@ARGV > 7 && $ARGV[7] eq "old") {
	$rna_seq_old = 1;
}
my @split_path_glassome = split('/', $path_glassome);
#Now add the last directory to the tscc path
$path_tscc .= "/" . $split_path_glassome[-1] . "/";

open OUT, ">tag_dir.sh";

print OUT "#!/bin/bash\n";
print OUT "#PBS -q hotel\n";
print OUT "#PBS -N " . $path_tscc . "\n";
print OUT "#PBS -l nodes=1:ppn=8\n";
print OUT "#PBS -l walltime=6:00:00\n";
print OUT "#PBS -o " . $path_tscc . "/output.txt\n";
print OUT "#PBS -e " . $path_tscc . "/error.txt\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M n\n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A glass-group\n";
print OUT "\n";
print OUT "cd " . $path_tscc . "\n";


#Create a processed folder
my $command = "mkdir " . $path_tscc . "/processed";
#print OUT $command . "\n";
`$command`;

$command = "mkdir " . $path_tscc . "/sam";
#print OUT $command . "\n";
`$command`;

#Move all sam files
$command = "mv " . $path_tscc . "*/*sam " . $path_tscc . "/processed";
`$command`;
#Create tag directories

my $output_dir = "";
my @a;
my @files = `ls $path_tscc/processed/*sam`;
my %atac = ();

foreach my $f (@files) {
	chomp $f;
	@a = split('\.', basename($f));
	$output_dir = $path_tscc . "/processed/" . $a[0];
	print OUT "makeTagDirectory " . $output_dir . " " . $f;
	if($method eq "STAR") {
		print OUT " -format sam";
		if($rna_seq_old == 0) {
			print OUT "  -flip";
		}
	} else {
		print OUT " -format sam";
		if(index(uc($f),"ATAC") != -1) {
			$atac{$a[0]} = 1;
		}

	}
	if($genome eq "mm10" || $genome eq "mm9" || $genome eq "hg19" || $genome eq "hg18") {
		print OUT " -genome $genome -checkGC\n";
	} else {
		print OUT "\n";
	}
}
use Data::Dumper;
#Copy all log files to tag directores
my @directory_name;
if($method eq "STAR") {
	@files = `ls $path_tscc/*/Log.final.out`;
	foreach my $f (@files) {
		chomp $f;
		@a = split('\.', basename($f));
	#	$output_dir = $path_tscc . "/processed/" . $a[0];
		@directory_name = split('\/', $f);
	#	print $f . "\n";
	#	print Dumper @directory_name;
		$output_dir = $directory_name[-2];
		print $output_dir . "\n";
		$command = "mv " . $f . " " . $path_tscc . "/processed/" . $output_dir . "/" . $output_dir . ".mapping.log";
		print OUT $command . "\n";
	}
} else {
	@files = `ls $path_tscc/*/*log`;
	foreach my $f (@files) {
		chomp $f;
		@a = split('\.', basename($f));
	#	$output_dir = $path_tscc . "/processed/" . $a[0];
		$output_dir = $a[0];
		$command = "mv " . $f . " " . $path_tscc . "/processed/" . $output_dir;
		print OUT $command . "\n";
	}
}
#Remove all sam files
$command = "mv " . $path_tscc . "/processed/*sam " . $path_tscc . "/sam/";
print OUT $command . "\n";

foreach my $files (keys %atac) {
	print OUT "rm " . $path_tscc . "/processed/" . $files . "/chrM.tags.tsv\n";
	print OUT "makeTagDirectory " . $path_tscc . "/processed/" . $files . "_without_M -d " . $path_tscc . "/processed/" . $files . "\n";
	print OUT "mv " . $path_tscc . "/processed/" . $files . "/tagInfo.txt " . $path_tscc . "/processed/" . $files . "_without_M/tagInfo_with_M.txt\n";
	print OUT "mv " . $path_tscc . "/processed/" . $files . "/*log " . $path_tscc . "/processed/" . $files . "_without_M/\n";
	print OUT "rm -rf " . $path_tscc . "/processed/" . $files . "\n";
}

#Move all tag directories to processed folder in archive
$command = "mkdir " . $path_glassome . "/processed";
print OUT $command . "\n";

$command = "cp -R " . $path_tscc . "/processed/Sample_* " . $path_glassome . "/processed";
print OUT $command . "\n";
close OUT;
`qsub tag_dir.sh`; 

