#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";
use Spiders::DatabaseQuery;
use Spiders::Params;
use Cwd 'abs_path';
use Data::Dumper;
use File::Basename;
use List::Util qw(min max sum);

my ($paramFile, $ms2OutFile) = @ARGV;
my $ms2File = $ms2OutFile;
$ms2File =~ s/\.out$//;

##########################################
## Parameter loading and initialization ##
##########################################
my %params = getParams($paramFile);
my $database = $params{'database'};
my $decoyDatabase = $database . "_decoy";
if ($database =~ /,/) {
	$decoyDatabase =~ s/,/_decoy,/;
}
my $neutralLossFile = $Bin . "/neutralLoss.csv";	
my $bondEnergyFile =  $Bin . "/bondenergies.txt";
my $H = 1.007276466812;

####################################
## Read .MS2 file to obtain peaks ##
####################################
open (MS2, "<", $ms2File) or die "Cannot open $ms2File\n";
my $header = <MS2>;
my (@mzArray, @intArray);
while (<MS2>) {
	chomp ($_);
	my ($mz, $intensity) = split(/\t/, $_);
	push (@mzArray, $mz);
	push (@intArray, $intensity);
}
close (MS2);

####################################
## Calculation of matching scores ##
####################################
## Choose a portion of high intensity peaks (measured/experimental/observed)
my $nPeaks = int(scalar(@mzArray) * ($params{'percentage_ms2_peaks'} / 100));
my @index = sort {$intArray[$b] <=> $intArray[$a]} (0..$#intArray);
@intArray = @intArray[@index[0..$nPeaks]];
@mzArray = @mzArray[@index[0..$nPeaks]];

## Make a hash which has a key of integerized m/z value for speed-up
my %mzHash;
for (my $i = 0; $i < scalar(@mzArray); $i++) {
	$mzHash{int($mzArray[$i])}{$mzArray[$i]} = $intArray[$i];
}

## Read .out file to calculate scores and write the result to .score file
my $hasError = 0;
open (OUT, "<", $ms2OutFile) or die "Cannot open $ms2OutFile\n";
$header = <OUT>;
chomp ($header);
my $scoreFile = $ms2File . ".score";
open (SCORE, ">", $scoreFile) or die "Cannot open $scoreFile\n";
$header = $header . "\t#measuredPeaks\t#theoreticalPeaks\t#matchedPeaks\tmscore\tIntensityRatio\n";
print SCORE $header;
while (<OUT>) {
	my $line = $_;
	chomp ($line);
	my ($index, $monoMass, $formula, $name, $struct, $inchi, $type, $adduct) = split(/\t/, $line);
	
	## Theoretical fragmentation
	
	
	## *************************************************** ##
	## Currently, only METFRAG-based method is implemented ##
	## To-do: CFM-ID based method needs to be implemented  ##
	## *************************************************** ##
	
	my @fragResult = smilesFragmentation($struct, 50, 2, "no", "yes", $neutralLossFile, $bondEnergyFile);
	if ($fragResult[0] =~ /Error/) {
		$hasError = 1;		
	}
	if ($hasError == 0) {
		my (@fragMzArray, @fragSmilesArray);
		for (my $i = 2; $i < scalar(@fragResult); $i++) {
			next if ($fragResult[$i] =~ /atom/);
			my @elems = split(/ /, $fragResult[$i]);
			next if($elems[3] !~ /\d+/);
			push (@fragMzArray, $elems[3]);
			push (@fragSmilesArray, $elems[0]);
		}
		
		## Manipulate the theoretically fragmented peaks according to the mode
		if ($params{'mode'} == 1) {	## Positive mode
			@fragMzArray = map {$_ - 0.00054858} @fragMzArray;
			if (lc($type) eq "decoy") {
				## Add one more H+ to all peaks (is it really a default decoy setting?)
				@fragMzArray = map {$_ + $H} @fragMzArray;
			}
		} elsif ($params{'mode'} == -1) {	## Negative mode
			@fragMzArray = map {$_ - 2 * $H + 0.00054858} @fragMzArray;
			if (lc($type) eq "decoy") {
				## Add one more H+ to all peaks (is it really a default decoy setting?)
				@fragMzArray = map {$_ + $H} @fragMzArray;
			}
		}
		if (defined ($adduct) && $adduct ne "NA") {
			$adduct = "adduct_" . $adduct;;
			@fragMzArray = map {$_ + $params{$adduct}} @fragMzArray;
		}
		
		## Calculate matching scores
		my ($nMatches, $matchedHash) = compareSpectra(\%mzHash, \@fragMzArray, \%params);
		my $pvalue = 1;
		my $intRatio = 0;
		my $mscore = 0;
		if ($nMatches > 0) {
			$pvalue = calcPvalue(\@mzArray, \@fragMzArray, $nMatches, \%params, $line);
			foreach my $matchedMz (keys %$matchedHash) {
				$intRatio += $mzHash{int($matchedMz)}{$matchedMz};
			}
			$intRatio /= sum(@intArray);
			$mscore = -log($pvalue) / log(10);
		} 
		print SCORE $line . "\t" . scalar(@mzArray) . "\t" . scalar(@fragMzArray) . "\t$nMatches\t$mscore\t$intRatio\n";
	} else {
		print SCORE $line . "\t" . scalar(@mzArray) . "\t0\t0\t0\t0\n";
	}
}
close (SCORE);
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
		if ($param =~/adduct\_(\w+)/) {
			$adducts{$1} = $$paramHash{$param};
		}
	}
	return \%adducts;
}

sub smilesFragmentation {
	my ($smiles, $mass, $depth, $breakAromaticRingFlag, $molecularFormulaRedundancyCheck, $neutralLossFile, $bondEnergyFile) = @_;
	if (!$smiles) {
		print "Formula is not defined\n";
		exit;
	}
	if (!$mass) {
		$mass = "0.0";
	}
	if (!$depth) {
		$depth = "2";
	}
	if (!$molecularFormulaRedundancyCheck) {
		$molecularFormulaRedundancyCheck = "no";
	}
	if (!$breakAromaticRingFlag) {
		$breakAromaticRingFlag = "yes";
	}	
	my $jarpath = $Bin . "/JumpJarPackage/JumpPackage.jar";
	$smiles =~ s/\(/\\\(/g;
	$smiles =~ s/\)/\\\)/g;
	$smiles ="\"" . $smiles . "\"";	
	my $command = "java -jar " . $jarpath . " -AIMFragmenterWrapper " . "\"java -jar " . $jarpath . "\" .tmp " . $smiles . " " . $mass . " " . $depth . " " . $breakAromaticRingFlag . " " . $molecularFormulaRedundancyCheck . " " . $neutralLossFile . " " . $bondEnergyFile;
    my @res = `$command`;
    return (@res);
}

sub compareSpectra {
	my ($measuredHash, $theoretical, $parameter) = @_;	
	my $tol = $$parameter{'frag_mass_tolerance'};
	my $tolUnit = 1;
	if (defined($$parameter{'frag_mass_tolerance_unit'})) {
		$tolUnit = $$parameter{'frag_mass_tolerance_unit'};
	}
	my %matchedHash;
	for (my $i = 0; $i < scalar(@$theoretical); $i++) {	
		my $theoMz = $$theoretical[$i];		
		if ($tolUnit == 2) {
			$tol = $tol * int($theoMz) / 1e6;
		}
		my $lL = int($theoMz - $tol);
		my $uL = int($theoMz + $tol) + 1;	
		for (my $intTheoMz = $lL; $intTheoMz < $uL; $intTheoMz++) {	## $intTheoMz = integerized theoretical m/z value
			if (defined($$measuredHash{$intTheoMz})) {
				foreach my $measuredMz (keys %{$$measuredHash{$intTheoMz}}) {
					if (($measuredMz + $tol) > $theoMz and ($measuredMz - $tol) < $theoMz) {
						next if ($matchedHash{$measuredMz});
						$matchedHash{$measuredMz} = 1;
					}
				}
			}
		}
	}
	my $nMatches = scalar (keys %matchedHash);
	return ($nMatches, \%matchedHash);
}

sub calcPvalue {
	my ($measuredMzArray, $fragMzArray, $nMatches, $params, $info) = @_;
	## Hypergeometric p-value is calculated using
	## n: # of all possible peaks (considering m/z tolerance)
	## k: # of measured peaks
	## r: # of theoretical peaks
	## x: # of matched peaks
		
	my $tol = $$params{'frag_mass_tolerance'};
	my $tolUnit = 1;
	if (defined($$params{'frag_mass_tolerance_unit'})) {
		$tolUnit = $$params{'frag_mass_tolerance_unit'};
	}	
	if ($tolUnit == 2) {
		## Use the average mass to convert the unit	
		$tol = $tol * ((min(@$measuredMzArray) + max(@$measuredMzArray)) / 2) / 1e6;
	}
	my $n = int( (max(@$measuredMzArray) - min(@$measuredMzArray)) / (2 * $tol) );
	my $k = scalar(@$measuredMzArray);
	my $r = scalar(@$fragMzArray);
	my $x = $nMatches;
	## Check: $x should be less than or equal to min($k, $r)
	my $pvalue = 1;
	if ($x <= min($r, $k)) {
		$pvalue = fisherExact($n, $k, $r, $x);
	} else {
		print "  Warning: unreasonable numbers of peaks. Please check the following entry\n";
		print "  $info\n";
	}
	return ($pvalue);
}

sub fisherExact {
	my ($n, $k, $r, $x) = @_;
	my @pr;
	for (my $i = 0; $i < $x; $i++) {
		if (($k - $i) >= 0 && ($r - $i) >= 0 && (($n - $r) - ($k - $i)) >= 0) {			
			push (@pr, hypergeometric($n, $k, $r, $i));
		} else {
			next;
		}
	}
	## One-tailed p-value
	my $p = 1 - sum(map {sprintf("%e", $_)} @pr);	## Scientific unit	
	if ($p <= 0) {
		$p = 1 - sum(map {sprintf("%.20f", $_)} @pr);	## 20 decimal points
	}
	if ($p <= 0) {
		$p = 1e-20;	## If p-value is still 0 or negative, then saturate to 20
	}	
	return ($p);
}

sub hypergeometric {
	my ($n, $k, $r, $x) = @_;
	my $pr = exp(choose($r, $x) + choose($n - $r, $k - $x) - choose($n, $k));
	return ($pr);
}

sub choose {
	my ($n, $k) = @_;
	my ($result, $j) = (0, 1);
	return 0 if $k > $n || $k < 0;
	$k = ($n - $k) if ($n - $k) < $k;
	while ($j <= $k) {
		$result += log($n--);
		$result -= log($j++);
	}
	return ($result);
}