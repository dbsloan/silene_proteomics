#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 input_peptide_file\n\n";

my $peptide_file = shift or die ($usage);

my @peptide_lines = file_to_array($peptide_file);


shift @peptide_lines;

#keep track of abundances by protein for each of the 6 samples, either without peaks classified as "Peak Found";
my %abundance_HoA;
my %abundance_noPF_HoA;

#keep track of unique peptides and total PSMs (across all samples) for each protein
my %unique_peptides;
my %total_psms;

foreach (@peptide_lines){
	my @sl = split (/\t/, $_);
	$sl[4] > 1 and next; #skip peptides that map to multiple proteins
	$sl[31] eq "Shared" and next; #this should be redundant with previous line but including it in case
	$sl[31] eq "No Quan Values" and next; #skip lines without reported values
	
	++$unique_peptides{$sl[7]};
	$total_psms{$sl[7]} += $sl[6];
			
	unless (exists$abundance_HoA{$sl[7]}){
		for (my $i=0; $i < 6; ++$i){
			$abundance_HoA{$sl[7]}[$i] = 0;
			$abundance_noPF_HoA{$sl[7]}[$i] = 0;
		}
	}
	
	for (my $i=0; $i < 6; ++$i){
		if ($sl[25+$i]){
			$abundance_HoA{$sl[7]}[$i] += $sl[25+$i];
			if ($sl[32+$i] eq "High"){
				$abundance_noPF_HoA{$sl[7]}[$i] += $sl[25+$i];			
			}
		}
	}
}

print "Protein\tUniquePeptides\tTotalPSMs\tCp1\tCp2\tLf1\tLf2\tMt1\tMt2\tCp1_noPF\tCp2_noPF\tLf1_noPF\tLf2_noPF\tMt1_noPF\tMt2_noPF\n";

foreach (sort keys %unique_peptides){
	
	print "$_\t$unique_peptides{$_}\t$total_psms{$_}";
	for (my $i=0; $i < 6; ++$i){
		print "\t$abundance_HoA{$_}[$i]";
	}
	for (my $i=0; $i < 6; ++$i){
		print "\t$abundance_noPF_HoA{$_}[$i]";
	}
	print "\n";
}

