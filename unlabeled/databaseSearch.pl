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

my $H = 1.007276466812;
my ($paramFile, $ms2File) = @ARGV;

##########################################
## Parameter loading and initialization ##
##########################################
my %params = getParams($paramFile);
my $database = $params{'database'};
my $decoyDatabase = $database . "_decoy";
if ($database =~ /,/) {
	$decoyDatabase =~ s/,/_decoy,/;
}

#######################################
## Database search using a .MS2 file ##
#######################################
my $query = new Spiders::DatabaseQuery();
$query -> setBin($Bin);
open (MS2, "<", $ms2File);
my $header = <MS2>;
close (MS2);
my $MH = (split(/\t/, $header))[0];	## Either (M+H)+ or (M-H)-
my $charge = (split(/\t/, $header))[1];

## Calculate 'neutral' monoisotopic mass to make a query
## Mass .MS2 (i.e. dta-format) file is (M+H)+ for positive mode or (M-H)- for negative mode, where M = neutral mass
my $neutralMass;
# = monoMass = ($mz - $H) * $charge + $H;
if ($params{'mode'} == -1) {
	$neutralMass = $MH + $H;
} else {
	$neutralMass = $MH - $H;
}

## DB search for formula
my ($allFormula, $targetFormula, $decoyFormula) = 
$query -> queryMassWithDatabaseFilter($neutralMass, 
                                      $params{'formula_mass_tolerance_searching'},
                                      $params{'mass_formula_database'},
                                      "yes", $database, $decoyDatabase);

## DB search for structures	
my (@targetStructures, @decoyStructures);
## 1. If there's neither target nor decoy, then search using adduct
## 2. Otherwise, go to structure database search
if (scalar(@$targetFormula) == 0 && scalar(@$decoyFormula) == 0 && $params{'adduct'} == 1) {
	## Adduct-based search
	my ($targetFormulaRef, $decoyFormulaRef, $targetStructRef, $decoyStructRef) = 
	adductBasedSearch($query, $neutralMass, \%params, $H);
	$targetFormula = $targetFormulaRef;
	$decoyFormula = $decoyFormulaRef;
	@targetStructures = @$targetStructRef;
	@decoyStructures = @$decoyStructRef;
} else {
	@targetStructures = searchStructures($query, $targetFormula, \%params, $H, 0, "NA");
	@decoyStructures = searchStructures($query, $decoyFormula, \%params, $H, 1, "NA");		
}

###########################################
## Print out the result to a file (.out) ##
###########################################
my $index = 0;
if (scalar(@targetStructures) == 0 && scalar(@decoyStructures) == 0) {
	## If there's neither targets nor decoys, just finish this script
	exit;
} else {
	my $outFile = $ms2File . ".out";
	open (OUT, ">", $outFile);
	print OUT "Index\tNeutralMass\tFormula\tName\tStructure\tInChi\tType\tAdduct\n";
	if (scalar(@targetStructures) > 0) {
		foreach my $targetStructure (@targetStructures) {
			chomp ($targetStructure);
			$index++;
			my @elems = split(/\t/, $targetStructure);
			my ($formula, $inchi, $struct, $name, $adduct) = @elems;
			## When there is only one element, '1' will be removed
			## e.g. C6H12O4N1S1 -> C6H12O4NS
			$formula =~ s/(\D+)1(\D+)/$1$2/g;
			$formula =~ s/(\D+)1$/$1/g;
			my $type = "Target";
			print OUT "$index\t$neutralMass\t$formula\t$name\t$struct\t$inchi\t$type\t$adduct\n";
		}
	}
	if (scalar(@decoyStructures) > 0) {
		foreach my $decoyStructure (@decoyStructures) {
			chomp ($decoyStructure);
			$index++;
			my @elems = split(/\t/, $decoyStructure);
			my ($formula, $inchi, $struct, $name, $adduct) = @elems;
			$formula =~ s/(\D+)1(\D+)/$1$2/g;
			$formula =~ s/(\D+)1$/$1/g;			
			my $type = "Decoy";
			print OUT "$index\t$neutralMass\t$formula\t$name\t$struct\t$inchi\t$type\t$adduct\n";
		}
	}
	close (OUT);
}

#################
## Subroutines ##
#################
sub adductBasedSearch {
	my ($query, $monoMass, $params, $H) = @_;
	my (@allTargetFormula, @allDecoyFormula, @allTargetStructures, @allDecoyStructures);
	my $adductHash = getAdducts($params);	
	foreach my $adductName (keys %$adductHash) {		
		my $monoMass = $monoMass - $$adductHash{$adductName};	## Monoisotopic mass without the adduct
		## DB search for formula
		my ($allFormula, $targetFormula, $decoyFormula) = 
		$query -> queryMassWithDatabaseFilter($monoMass, 
		                                      $params{'formula_mass_tolerance_searching'},
		                                      $params{'mass_formula_database'},
		                                      "yes", $database, $decoyDatabase);
		my @targetStructures = searchStructures($query, $targetFormula, $params, $H, 0, $adductName);
		my @decoyStructures = searchStructures($query, $decoyFormula, $params, $H, 1, $adductName);
		push (@allTargetFormula, @$targetFormula);
		push (@allDecoyFormula, @$decoyFormula);
		push (@allTargetStructures, @targetStructures);
		push (@allDecoyStructures, @decoyStructures);
	}
	return (\@allTargetFormula, \@allDecoyFormula, \@allTargetStructures, \@allDecoyStructures);	
}

sub searchStructures {
	my ($query, $formulaArray, $params, $H, $isDecoy, $adduct) = @_;
	my $decoyH = 1;	## Number of hydrogen in decoy
	if (defined($$params{'decoy_strategy'})) {
		$decoyH = $$params{'decoy_strategy'};
	}
	my @structures;
	foreach my $elem (@$formulaArray) {
		chomp ($elem);
		my ($formula, $mass) = split(/:/, $elem);
		$formula =~ s/(\D+)1(\D+)/$1$2/g;
		$formula =~ s/(\D+)1$/$1/g;
		my ($queryFormula, $queryMass);
		if ($isDecoy == 1) {
			if ($formula =~ /(.*H)(\d+)(\D+.*)/) {
				my $Hminus = $2 - $decoyH;
				$queryFormula = $1 . $Hminus . $3;
			}
			$queryMass = $mass - $H * $decoyH;
		} else {
			$queryFormula = $formula;
			$queryMass = $mass;
		}
		$query -> setDatabase($$params{'structure_database'});
		$query -> setFormula($queryFormula);
		$query -> setMass($queryMass);
		$query -> setDatatype("INCHIKEY,SMILES,GENERALNAME");
		$query -> setDBname($$params{'database'});				
		my $structure = $query -> queryStructureDatabase();
		chomp (@$structure);	## Each line in @$structure has a new line, so it needs to be removed
		if ($isDecoy == 1) {
			## Although 'decoy' is being searched using 'queryFormula',
			## the original formula should be written in the output file
			for (my $i = 0; $i < scalar(@$structure); $i++) {
				my @elems = split(/\t/, $$structure[$i]);
				my $newStructure = join("\t", ($formula, @elems[1..$#elems]));
				$$structure[$i] = $newStructure;
			}
		}
		@$structure = map {$_ . "\t$adduct"} @$structure;
		push (@structures, @$structure);		
	}
	return (@structures);
}

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
	return \%adducts;
}