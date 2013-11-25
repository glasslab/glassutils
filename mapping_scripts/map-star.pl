#!/usr/bin/perl -w

use strict;
#Script to map with STAR

my $starDir = "/projects/glass-group/bioinformatics/STAR/";

my $genome = $ARGV[0];
my $file = $ARGV[1];

my @split = split('/', $file);

my $wd = $split[0];

for(my $i = 1; $i < @split - 2; $i++) {
	$wd = $wd . "/" . $split[$i]; 
}
my $command = $starDir . "STAR --genomeDir " . $starDir . "genomes/" . $genome . " --readFilesIn " . $file . " --runThreadN 4";
`$command`;

$command = "mv " . $wd . "/" . $split[-2] . "/Aligned.out.sam " . $wd . "/" . $split[-2] . "/" . $split[-1] . "." . $genome . ".STAR.sam";
print $command . "\n";
`$command`;
