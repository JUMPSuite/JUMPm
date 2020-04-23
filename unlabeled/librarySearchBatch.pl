#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin);
use List::Util qw(min max);
use List::MoreUtils qw(indexes first_index);
use Statistics::Lite qw(mean median);

##########################################
## Parameter loading and initialization ##
##########################################
#my $paramFile = @ARGV;
my $paramFile = "../IROA_HILIC_NEG_Target/jumpm_negative.params";
my %params = getParams($paramFile);
my $H = 1.007276466812;
my $matchMzTol = 10;  ## Unit of ppm
my $matchRtTol = 10;  ## Unit of second

############################
### Read library entities ##
############################
### Column names are determined based on 
### type of LC column (C4, C8, C18 or HILIC) and mode (pos or neg),
#print "Loading a library\n";
#my $columnInfo = lc($params{'LC_column'});
#if ($params{'mode'} == -1) {
#	$columnInfo .= "n";
#} elsif ($params{'mode'} == 1) {
#	$columnInfo .= "p";
#}
#my @libInfo;
#my $libFile = $params{'library'};
#open (LIB, "<", $libFile) or die "Cannot open $libFile\n";
#my $header = <LIB>;
#my @headerElems = split(/\t/, $header);
#my $i = 0;
#while (<LIB>) {
#	chomp($_);
#	my @elems = split(/\t/, $_);
#	@{$libInfo[$i]}{@headerElems} = @elems;
#
#	## Open .MS2 file and calculate m/z of the library feature
#	my $ms2File = $libInfo[$i]{"$columnInfo" . "_linkms2"};
#	open (MS2, "<", $ms2File) or die "Cannot open $ms2File\n";
#	my $header = <MS2>;
#	while (<MS2>) {
#		chomp ($_);
#		my ($mz, $int) = split(/\t/, $_);
#		push (@{$libInfo[$i]{"$columnInfo" . "_ms2"}{'mz'}}, $mz);
#		push (@{$libInfo[$i]{"$columnInfo" . "_ms2"}{'int'}}, $int);
#	}
#	close (MS2);
#	chomp($header);
#	my ($MH, $charge) = split(/\t/, $header);	## Note that .MS2 file has M+H (or M-H) 
#
#	## Check the charge state
#	if ($libInfo[$i]{"$columnInfo" . "_charge"} > 0 && $libInfo[$i]{"$columnInfo" . "_charge"} != $charge) {
#		die "Charge state of the $i-th library entry\n";
#	}
#
#	my $mz;
#	if ($params{'mode'} == 1) {
#		$mz =(($MH - $H) + $charge * $H) / $charge;
#	} elsif ($params{'mode'} == -1) {
#		$mz = (($MH + $H) - $charge * $H) / $charge;
#	}
#	$libInfo[$i]{"$columnInfo" . "_mz"} = $mz;
#	$i++;
#	print "\rParsing #$i entries in the library";
#}
#close (LIB);
#print "\n";

###########################
## Read library entities ##
###########################
## Column names are determined based on 
## type of LC column (C4, C8, C18 or HILIC) and mode (pos or neg),
print "Loading a library\n";


$params{'LC_column'} = "hilic";


my $columnInfo = lc($params{'LC_column'});
if ($params{'mode'} == -1) {
	$columnInfo .= "n";
} elsif ($params{'mode'} == 1) {
	$columnInfo .= "p";
}
my @libInfo;
my $libFile = "./IROAlibrary/Metabolome_library_v0.1.4.txt";
my $libMS2Path = "./IROAlibrary/MS2/";
open (LIB, "<", $libFile) or die "Cannot open $libFile\n";
my $header = <LIB>;
my @headerElems = split(/\t/, $header);
my $i = 0;
while (<LIB>) {
	chomp($_);
	my @elems = split(/\t/, $_);
	@{$libInfo[$i]}{@headerElems} = @elems;

	## Open .MS2 file and calculate m/z of the library feature
	my $ms2File = $libMS2Path . basename($libInfo[$i]{"$columnInfo" . "_linkms2"});
	open (MS2, "<", $ms2File) or die "Cannot open $ms2File\n";
	my $header = <MS2>;
	while (<MS2>) {
		chomp ($_);
		my ($mz, $int) = split(/\t/, $_);
		push (@{$libInfo[$i]{"$columnInfo" . "_ms2"}{'mz'}}, $mz);
		push (@{$libInfo[$i]{"$columnInfo" . "_ms2"}{'int'}}, $int);
	}
	close (MS2);
	chomp($header);
	my ($MH, $charge) = split(/\t/, $header);	## Note that .MS2 file has M+H (or M-H) 

	## Check the charge state
	if ($libInfo[$i]{"$columnInfo" . "_charge"} > 0 && $libInfo[$i]{"$columnInfo" . "_charge"} != $charge) {
		die "Charge state of the $i-th library entry\n";
	}

	my $mz;
	if ($params{'mode'} == 1) {
		$mz =(($MH - $H) + $charge * $H) / $charge;
	} elsif ($params{'mode'} == -1) {
		$mz = (($MH + $H) - $charge * $H) / $charge;
	}
	$libInfo[$i]{"$columnInfo" . "_mz"} = $mz;
	$i++;
	print "#$i\tmz = $mz\tcharge = $charge\n";
}
close (LIB);
print "\n";

##############################
## Read feature information ##
##############################
## Note that the header of .MS2 file contains (M+H)+ or (M-H)-, and charge state
print "\nLoading feature information\n";
my @featInfo;
my $featFile = "/Research/Projects/7Metabolomics/IROA_HILIC_NEG_Target/test_fully_aligned.feature"; 
my $ms2Path = "/Research/Projects/7Metabolomics/IROA_HILIC_NEG_Target/MS2/";
open (FEATURE, "<", $featFile) or die "Cannot open $featFile\n";
$header = <FEATURE>;
@headerElems = split(/\t/, $header);
my @colRT = indexes {$_ =~ /_RT/} @headerElems;
my @colIntensity = indexes {$_ =~ /_Intensity/} @headerElems;
my @colZ = indexes {$_ =~ /_z$/} @headerElems;
$i = 1;
while (<FEATURE>) {
	print "\rParsing #$i features";
	chomp ($_);
	my @elems = split(/\t/, $_);
	@{$featInfo[$i]}{@headerElems} = @elems;
	my $meanRt = mean(grep {$_ ne "NA"} @elems[@colRT]);	## Some columns have "NA". They should be removed
	my $meanIntensity = mean(grep {$_ ne "NA"} @elems[@colIntensity]);
	$featInfo[$i]{'meanRt'} = $meanRt;
	$featInfo[$i]{'meanIntensity'} = $meanIntensity;
	
	## Open .MS2 file and obtain the charge state
	my $ms2File = $ms2Path . "/f" . $i . "\.MS2";
	if (-e $ms2File) {
		open (MS2, "<", $ms2File) or die "Cannot open $ms2File\n";
		my $header = <MS2>;
		while (<MS2>) {
			chomp ($_);
			my ($mz, $int) = split(/\t/, $_);
			push (@{$featInfo[$i]{"ms2"}{"mz"}}, $mz);
			push (@{$featInfo[$i]{"ms2"}{"int"}}, $int);
		}	
		close (MS2);
		chomp ($header);
		my ($MH, $charge) = split(/\t/, $header);	## Note that .MS2 file has M+H (or M-H) 
		$featInfo[$i]{'charge'} = $charge;
		
		## Check the consistency of m/z
		my $mz;
		if ($params{'mode'} == 1) {
			$mz =(($MH - $H) + $charge * $H) / $charge;
		} elsif ($params{'mode'} == -1) {
			$mz = (($MH + $H) - $charge * $H) / $charge;
		}
		my $diff = ($mz - $featInfo[$i]{'meanMz'}) / $mz * 1e6;
		if ($diff > 20) {
			print "Too much different\n";
		}
	}
	$i++;	
}
close (FEATURE);

##################
## RT alignment ##
##################
## Create temporary files; one for RT of library metabolites and the other for RT of features
## Just based on m/z similarity, select putative matched features and then perform RT alignment
open (REF, ">", "refRt.txt");
open (COMP, ">", "compRt.txt");
for (my $i = 0; $i < scalar(@libInfo); $i++) {
	my $refMz = $libInfo[$i]{"$columnInfo" . "_mz"};
	my $refRt = $libInfo[$i]{"$columnInfo" . "_rt"} * 60;	## Convert minute to second
	my $refZ = $libInfo[$i]{"$columnInfo" . "_charge"};
	my ($compRt, $compInd);	
	my $compMaxInt = 0;
	for (my $j = 1; $j < scalar(@featInfo); $j++) {		
		my $compMz = $featInfo[$j]{'meanMz'};
		my $compInt = $featInfo[$j]{'meanIntensity'};
		my $compZ = $featInfo[$j]{'charge'};
		my $mzDiff = abs($refMz - $compMz) / $refMz * 1e6;
		if ($refZ == 0) {
#			if ($mzDiff < $matchMzTol && defined($compZ) && $compInt >= $compMaxInt) {
			if ($mzDiff < $matchMzTol && $compInt >= $compMaxInt) {
				$compRt = $featInfo[$j]{'meanRt'};
				$compInd = $j;
				$compMaxInt = $compInt;
			}
		} else {
#			if ($mzDiff < $matchMzTol && defined($compZ) && $refZ == $compZ && $compInt >= $compMaxInt) {
			if ($mzDiff < $matchMzTol && $compInt >= $compMaxInt) {
				$compRt = $featInfo[$j]{'meanRt'};
				$compInd = $j;
				$compMaxInt = $compInt;
			}
		}
	}
	if (defined ($compRt)) {
		print REF "$refRt\t$i\n";
		print COMP "$compRt\t$compInd\n";
	}
}
close (REF);
close (COMP);
open (COMP_NEW, ">", "compRt_new.txt");
for (my $i = 1; $i < scalar(@featInfo); $i++) {		
	my $compRt = $featInfo[$i]{'meanRt'};
	print COMP_NEW "$i\t$compRt\n";
}
close (COMP_NEW);

## Alignment operation in R
my $command = "Rscript $Bin/R/simpleAlignment.r refRt.txt compRt.txt compRt_new.txt";
system ($command);

## Update RTs to aligned ones
open (RES, "<", "alignment_result.txt");
while (<RES>) {
	chomp ($_);
	my ($ind ,$rtNew) = split(/\t/, $_);
	$featInfo[$ind]{'meanRtNew'} = $rtNew;
}
close (RES);

## Remove temporary files
$command = "rm refRt.txt compRt.txt compRt_new.txt alignment_result.txt";
system ($command);	
print "\n";

############################################
## Match features and library metabolites ##
############################################
my %res;
my $minP = 1;
open (ERR, ">", "measured_rt_error.txt");
for (my $i = 1; $i < scalar(@featInfo); $i++) {
	## Compare each feature's m/z and RT with every library metabolite
	for (my $j = 0; $j < scalar(@libInfo); $j++) {		
		## 1. Charge states should be identical(?)
		my $libZ = $libInfo[$j]{"$columnInfo" . "_charge"};
		my $featZ = $featInfo[$i]{'charge'};
		next if (!defined ($featZ));
		next if ($featZ > 0 && $libZ > 0 && $libZ != $featZ);
		## 2. Check m/z- and RT-differences
		my $featMz = $featInfo[$i]{'meanMz'};
		my $featRt = $featInfo[$i]{'meanRtNew'};
		my $libMz = $libInfo[$j]{"$columnInfo" . "_mz"};
		my $libRt = $libInfo[$j]{"$columnInfo" . "_rt"} * 60;	## Note that library RT is using "minute"
		my $mzDiff = abs($libMz - $featMz) / $libMz * 1e6;
		my $rtDiff = abs($featRt - $libRt);
		if ($mzDiff < $matchMzTol) {
			## Calculate the similarity between feature and library MS2 spectra
			print ERR "$i\t$j\t$rtDiff\n";
			my $sim = calcMS2Similarity($featInfo[$i]{"ms2"}, $libInfo[$j]{"$columnInfo" . "_ms2"});
			$res{$i}{$j}{'sim'} = $sim;
			$res{$i}{$j}{'simP'} = 1 - $sim;
			if ((1 - $sim) > 0) {
				$minP = min(1 - $sim, $minP);
			}
			$res{$i}{$j}{'rtErr'} = $rtDiff;
		}
	} 
} 
close (ERR);

## Calculate p-value-like value for RT-differences
$command = "Rscript $Bin/R/calcAlignmentPvalue.R alignment_residual.txt measured_rt_error.txt";
system ($command);

## Organize the output
open (RES, "<", "calculated_rt_pvalue.txt") or die "Cannot open calculated_rt_pvalue.txt";
while (<RES>) {
	chomp;
	my ($i, $j, $p) = split(/\t/, $_);
	$res{$i}{$j}{'rtP'} = $p;
	if ($p > 0) {
		$minP = min($minP, $p);
	}
}
close (RES);

open (OUTPUT, ">", "libSearchResult.txt");
print OUTPUT "No\tFeature_Index\tFeature_m\/z\tFeature_RT(original)\tFeature_RT(aligned)\tFeature_Intensity\t" . 
			"Formula\tName\tSMILES\tInChiKey\tRT_shift\tMS2_similarity\tScore\n";
my $index = 1;
foreach my $key1 (keys %res) {
	foreach my $key2 (keys %{$res{$key1}}) {
		## key1 = feature number/index
		## key2 = library entity index		
		my $featMz = $featInfo[$key1]{'meanMz'};
		my $featOrigRt = $featInfo[$key1]{'meanRt'};
		my $featAlignedRt = $featInfo[$key1]{'meanRtNew'};
		my $featIntensity = $featInfo[$key1]{'meanIntensity'};
		my $libFormula = $libInfo[$key2]{'formula'};
		my $libName = $libInfo[$key2]{'name'};
		my $libSmiles = $libInfo[$key2]{'SMILES'};
		my $libInChI = $libInfo[$key2]{'InChIKey'};
		my $rtErr = $res{$key1}{$key2}{'rtErr'};
		my $rtP = $res{$key1}{$key2}{'rtP'};
		if ($rtP == 0) {
			$rtP = $minP / 10;
		}
		my $sim = $res{$key1}{$key2}{'sim'};
		my $simP = $res{$key1}{$key2}{'simP'};
		if ($simP == 0) {
			$simP = $minP / 10;
		}
		my $score = -2 * (log($rtP) + log($simP));	## Fisher's method for combining p-values
		print OUTPUT "$index\t$key1\t$featMz\t$featOrigRt\t$featAlignedRt\t$featIntensity\t$libFormula\t$libName\t$libSmiles\t$libInChI\t$rtErr\t$sim\t$score\n";
		$index++;
	}
}
close (OUTPUT);

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

sub calcMS2Similarity {
	## Clustering millions of tandem mass spectra, J Proteome Res. 2008; 7: 113-22
	my ($feat, $lib) = @_;	
	my $k = min(scalar(@{$$feat{'mz'}}), scalar(@{$$lib{'mz'}}));
	$k = min($k, 30);
	
	## Sort arrays in descending order of intensity and keep $k strongest peaks
	my @ind = sort {$$feat{'int'}[$b] <=> $$feat{'int'}[$a]} 0..$#{$$feat{'int'}};
	@{$$feat{'mz'}} = @{$$feat{'mz'}}[@ind[0..($k - 1)]];
	@{$$feat{'int'}} = @{$$feat{'int'}}[@ind[0..($k - 1)]];
	@ind = sort {$$lib{'int'}[$b] <=> $$lib{'int'}[$a]} 0..$#{$$lib{'int'}};
	@{$$lib{'mz'}} = @{$$lib{'mz'}}[@ind[0..($k - 1)]];
	@{$$lib{'int'}} = @{$$lib{'int'}}[@ind[0..($k - 1)]];
	my (%featHash, %libHash);
	for (my $i = 0; $i < $k; $i++) {
		$featHash{$$feat{'mz'}[$i]} = $$feat{'int'}[$i];
		$libHash{$$lib{'mz'}[$i]} = $$lib{'int'}[$i];
	}
	
	## Join two sets of m/z and make a new set of unique m/z
	## Duplicate masses are removed as follows
	## - We consider two peaks to have a similar mass if they are within 0.5 Da from each other)
	## - For those two peaks having similar mass, the lower one will be the unique one 
	##   (e.g. One peak with m/z = 100 and the other peak with m/z = 100.4 -> they will be merged to m/z = 100)
	my @mzArray = (@{$$feat{'mz'}}, @{$$lib{'mz'}});
	@mzArray = sort {$a <=> $b} @mzArray;
	my %mzHash;
	my $val = 0;
	foreach my $mz (@mzArray) {
		if (abs($mz - $val) <= 0.5) {
			$mzHash{$mz} = $val;
		} else {
			$mzHash{$mz} = $mz;
			$val = $mz;
		}
	}
	
	## Reduction of spectrum to a vector by assigning to each intensity to the unique m/z bins
	## And then, calculate the similarity; normalized dot-product
	## Initialization
	my %s;
	while (my ($key, $val) = each %mzHash) {
		$s{$val}{'feat'} = 0;
		$s{$val}{'lib'} = 0;
	}
	
	while (my ($key, $val) = each %mzHash) {
		if (defined $featHash{$key}) {
			$s{$val}{'feat'} += sqrt($featHash{$key}); 
		}
		if (defined $libHash{$key}) {
			$s{$val}{'lib'} += sqrt($libHash{$key});
		}
	}
	my ($den, $num1, $num2) = (0, 0, 0);
	foreach my $mz (keys %s) {
		$den += $s{$mz}{'feat'} * $s{$mz}{'lib'};
		$num1 += $s{$mz}{'feat'} ** 2;
		$num2 += $s{$mz}{'lib'} ** 2;
	}
	my $sim = $den / sqrt($num1 * $num2);
	return ($sim);
}

