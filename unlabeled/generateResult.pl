#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

## Initialization
my ($ms2Path, $featureFile, $paramFile) = @ARGV;
#my $ms2Path = "/Research/Projects/7Metabolomics/JUMPm/IROAsamples/20191126/align_IROA_IS_NEG.1";
#my $featureFile = "./IROAsamples/20191126/IROA_IS_NEG_fully_aligned.feature";
#my $paramFile = "./IROAsamples/20191126/jumpm_negative.params";
my $H = 1.007276466812;

## Parse parameters/adduct information
my %params = getParams($paramFile);
my %adducts = getAdducts(\%params);

## Loading features
my %featureHash;
open (FEATURE, "<", $featureFile) or die "Cannot open $featureFile\n";
my $featureHeader = <FEATURE>;
my (@tmp, @intColInd);
my $colInd = 0;
foreach (split("\t", $featureHeader)) {
	if (lc($_) =~ /intensity/) {
		push (@tmp, $_);
		push (@intColInd, $colInd);
	}
	$colInd++;
}
$featureHeader = join("\t", @tmp);
while (<FEATURE>) {
	chomp($_);
	my @elems = split("\t", $_);
	my ($num, $ion, $mz, $rt) = @elems[0..3];
	$featureHash{$num}{'ion'} = $ion;
	($featureHash{$num}{'charge'}) = $ion =~ /\[*.\](\d)[+-]/;
	if (!defined $featureHash{$num}{'charge'}) {
		$featureHash{$num}{'charge'} = 1;
	}
	$featureHash{$num}{'mz'} = $mz;
	$featureHash{$num}{'rt'} = $rt;
	$featureHash{$num}{'intensity'} = join("\t", @elems[@intColInd]);
}
close (FEATURE);

## Loading .score files
$ms2Path =~ s/\/$//;
my $outFile = $ms2Path . ".result";
my $header = "FeatureNo\tIon\tm\/z\tRT\tFormula\tTarget\/Decoy\t" . 
			"Name\tSMILES\tInChiKey\tMscore\t" . $featureHeader;
open (OUT, ">", $outFile) or die "Cannot open $outFile\n";
print OUT $header;
## Sort .score files according to the feature number
my @scoreFiles = glob("$ms2Path/*.score");
my @featureNums;
foreach my $scoreFile (@scoreFiles) {
	my ($featureNum) = $scoreFile =~ /f(\d+)\.MS2\.score/;
	push (@featureNums, $featureNum); 
}
my @index = sort {$featureNums[$a] <=> $featureNums[$b]} (0..$#featureNums);
@featureNums = @featureNums[@index];
@scoreFiles = @scoreFiles[@index];

## Read each .score file and write to output files
for (my $i = 0; $i < scalar(@scoreFiles); $i++) {
	open (SCORE, "<", $scoreFiles[$i]) or die "Cannot open $scoreFiles[$i]\n";	
	my ($featureNo) = basename($scoreFiles[$i]) =~ /f(\d+).MS2.score/;
	my (@entries, @mscores);
	while (<SCORE>) {
		chomp ($_);		
		next if ($_ =~ "Index");
		my @elems = split(/\t/, $_);
		my ($dummy1, $neutralMass, $formula, $name, $inchi, $inchikey, $type, $adduct, $dummy2, $dummy3, $dummy4, $mscore, ) = @elems;
		my $ion = $featureHash{$featureNo}{'ion'};
		my $charge = $featureHash{$featureNo}{'charge'};
		## Handling of adducts
		if ($adduct ne "NA") {
			if ($adduct eq "2H") {
				my $coeff = $charge + 1;
				$ion =~ s/(M-.*)H/M-\Q$coeff\EH/;				
			} elsif ($adduct eq "-2H+Na") {
				my $coeff = $charge + 1;
				$ion =~ s/(M-.*)H/M-\Q$coeff\EH/;
				$ion =~ s/(M[+-].*H)/$1+Na/;
			} elsif ($adduct eq "-2H+K") {
				my $coeff = $charge + 1;
				$ion =~ s/(M-.*)H/M-\Q$coeff\EH/;
				$ion =~ s/(M[+-].*H)/$1+K/;
			} else {
				$ion =~ s/(M[+-].*H)/$1+$adduct/;
			}
		}
		my $mz = $featureHash{$featureNo}{'mz'};
		my $rt = $featureHash{$featureNo}{'rt'};
		my $intensity = $featureHash{$featureNo}{'intensity'};
		my $line = "$featureNo\t$ion\t$mz\t$rt\t$formula\t$type\t$name\t$inchi\t$inchikey\t$mscore\t$intensity"; 
		push (@entries, $line);
		push (@mscores, $mscore);
	}
	close (SCORE);
	my @ind = sort {$mscores[$b] <=> $mscores[$a]} (0..$#mscores);
	@entries = @entries[@ind];
	for (my $j = 0; $j < scalar(@entries); $j++) {
		print OUT "$entries[$j]\n";
	}
}
close (OUT);

#################
## Subroutines ##
#################
sub getParams {
	my ($file) = @_;
	if (!-e $file) {
		die "Cannot open a parameter file. Please check the path of the parameter file\n";
	}
	my %params;
	my $nInputFiles = 0;
	open (PARAMS, "<", $file);
	while (<PARAMS>) {
		next if (!/\=/ || /^\#/ || /^\s+\#/);
		$_ =~ s/\s+//g;
		$_ =~ s/\#.*//;
		my ($key, $value)  = split("\=", $_);
		if ($key eq "input_file") {
			$params{$key}[$nInputFiles] = $value;
			$nInputFiles++;
		} else {
			$params{$key} = $value;
		}
				
	}
	close (PARAMS);
	return (%params);
}

sub getAdducts {
	my ($paramHash) = @_;
	my %adducts;
	foreach my $param (keys %$paramHash) {
		if ($param =~/adduct\_(.*)/) {
			$adducts{$1} = $$paramHash{$param};
		}
	}
	return (%adducts);
}