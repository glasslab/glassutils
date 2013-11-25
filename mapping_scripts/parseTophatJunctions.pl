#!/usr/bin/perl -w
if (@ARGV < 1) {
	print STDERR "<tophat junctions.bed file>\n";
	exit;
}
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^track/);
	my @line = split /\t/;

	my @offset = split /\,/,$line[10];
	my $s = $line[1]+$offset[0];
	my $e = $line[2]-$offset[1]+1;
	print "$line[0]\t$s\t$e\t$line[5]\n";
}
