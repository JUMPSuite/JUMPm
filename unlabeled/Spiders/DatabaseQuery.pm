#!/usr/local/bin/perl

package Spiders::DatabaseQuery;
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";
use POSIX;
no warnings 'uninitialized';
use vars qw($VERSION @ISA @EXPORT);

$VERSION = 1.0;
@ISA = qw(Exporter);
@EXPORT = qw(setBin getBin setDatabase getDatabase setFormula getFormula setDatatype getDatatype setDBname getDBname setMass getMass setParameter getParameter setDtaFile getDtaFile queryMassWithDatabase queryMassWithDatabaseFilter queryStructure formulaParser queryStructureDatabase AIMPuritySampler  partialCandidate getJarPath);
 
sub new {
	my ($class, %arg) = @_;
	my $self = {};
	bless $self, $class;
	return $self;
}

sub setBin {
	my ($self, $Bin)=@_;
	$self -> {_Bin} = $Bin;
}

sub getBin {
	my ($self) = @_;
	return $self -> {_Bin};	
}

sub setDatabase {
	my ($self, $database) = @_;
	$self -> {_database} = $database;	
}

sub getDatabase {
	my ($self) = @_;
	return $self -> {_database};	
}

sub setFormula {
	my ($self, $formula) = @_;
	$self -> {_formula} = $formula;	
}

sub getFormula {
	my ($self) = @_;
	return $self -> {_formula};	
}

sub setDatatype {
	my ($self, $datatype) = @_;
	$self -> {_datatype} = $datatype;	
}

sub getDatatype {
	my ($self) = @_;
	return $self -> {_datatype};	
}

sub setDBname {
	my ($self, $DBname) = @_;
	$self -> {_dbname} = $DBname;	
}

sub getDBname {
	my ($self) = @_;
	return $self -> {_dbname};	
}

sub setMass {
	my ($self, $mass) = @_;
	$self -> {_mass} = $mass;	
}

sub getMass {
	my ($self) = @_;
	return $self -> {_mass};	
}

sub setParameter {
	my ($self, $param) = @_;
	$self -> {'_parameter'} = $param;
}

sub getParameter {
	my $self = shift;
	return $self -> {'_parameter'};
}

sub setDtaFile {
	my ($self, $dtaFile) = @_;
	return $self -> {_dta_file} = $dtaFile;
}

sub getDtaFile {
	my $self = shift;
	return $self -> {_dta_file};}


#########################################
## Query the formula with a given mass ##
#########################################
sub queryMassWithDatabase {
	my ($self, $mass, $ppm, $isotopeDatabase, $rgdbCheck, $dbName) = @_;
	if (!$ppm) {
		$ppm = 1;
	}
	if (!$mass) {
		print "Mass is not defined\n";
		return ();
	}
	if (!$isotopeDatabase) {
		print "IsotopeDatabase path is not defined\n";
		return ();	
	}
	if (!$rgdbCheck) {
		$rgdbCheck = "no";
	}
	my $jarPath = $self -> getJarPath();
	my $command = "java -jar " . $jarPath . " -queryMassWithPPM " . $mass . " " . 
				   $ppm . " " . $isotopeDatabase . " " . $rgdbCheck . " " . $dbName;
    my @value = `$command`;
	return \@value;
}

sub queryMassWithDatabaseFilter {
	my ($self, $mass, $ppm, $isotopeDatabase, $rgdbCheck, $targetDB, $decoyDB) = @_;
	if (!$ppm) {
		$ppm = 1;
	}
	if (!$mass) {
		print "mass is not defined\n";
		return ();
	}
	if (!$isotopeDatabase) {
		print "IsotopeDatabase path is not defined\n";
		return ();   
	}
	if (!$rgdbCheck) {
		$rgdbCheck = "no";
	}
	if (!$targetDB) {
		print "TargetDB term is not defined\n";
		return ();
	}
	if (!$decoyDB) {
		print "DecoyDB term is not defined\n";
		return ();
	}
	my $jarPath = $self -> getJarPath();
	my $command = "java -jar " . $jarPath . " -QueryMassWithPPMFilter " . $mass . " " . 
				   $ppm . " " . $isotopeDatabase . " " . $rgdbCheck . " " . $targetDB . " " . $decoyDB;
	my @value = `$command`;
	
	## Organize the search result
	my (%all, %target, %decoy);
	foreach (@value) {
		my $line = $_;
		my @tags = split('\t', $line);
		if ($tags[0] eq "All") {
			$all{$tags[1]} = $tags[1];
		} elsif ($tags[0] eq "Target") {
			$target{$tags[1]} = $tags[1];
		} elsif ($tags[0] eq "Decoy") {
			$decoy{$tags[1]} = $tags[1];
		} else {
			print "This shouldn't happen\n";
		}
	}
	my @allArray;
	for (keys %all) {
		push (@allArray, $_);
	}
	my @targetArray;
	for (keys %target) {
		push (@targetArray, $_);
	}
	my @decoyArray;
	for (keys %decoy) {
		push (@decoyArray, $_);
	}
	return (\@allArray, \@targetArray, \@decoyArray);
}

sub queryStructure {
	my ($self) = @_;
	my @str = ();
	my $database = $self -> getDatabase();
	my $mass = $self -> getMass();
	my $formula = $self -> getFormula();
	my $formulaFile = $formula;
	$formulaFile =~ s/\d+//g;
	my ($massInteger, $massDecimal) = split(/\./, $mass);
	my $structureFile = "$database/$massInteger/$formulaFile" . ".txt";
	if (-e($structureFile)) {
		open (FILE, "<", $structureFile);
		while (<FILE>) {
			my @data = split(/\t/, $_);
			next if (!defined($formula));
			next if (!defined($data[1]));
			my $element1 = $self -> formulaParser($formula);
			my $element2 = $self -> formulaParser($data[1]);			
			my $formula1 = "";
			my $formula2 = "";
			foreach my $key (sort keys %$element1) {
				my $elm = $key . $element1 -> {$key};
				$formula1 .= $elm;
			}
			foreach my $key (sort keys %$element2) {
				my $elm = $key . $element2 -> {$key};
				$formula2 .= $elm;
			}			
			if ($formula1 eq $formula2) {			
				my $smiles = join("\t", $data[2], @data[4..5]);
				push (@str, $smiles);
			}
		}
	}
	return \@str;
}

sub formulaParser {
	my ($self, $formula) = @_;
	my %elements;
	while ($formula =~ m/([A-Z][a-z]*)(\d+\.?\d*)?/go) { # XXX new
		my ($elm, $count) = ($1, $2);
		$count = 1 unless defined $count;
		if (defined $elements{$elm}) {
			$elements{$elm} += $count;
		} else {
			$elements{$elm} = $count;
		}
	}
    return \%elements;
}

#################################################################################
## Query the structure given a formula                                         ##
## $datatype = INCHIKEY;INCHISTR;SMILE,IUPAC,GENERALNAME, (multiple selection) ##
## $dbname = PUBCHEM,HMDB,KEGG                                                 ##
#################################################################################
sub queryStructureDatabase {
	my $self = shift;
	my $database = $self -> getDatabase();
	my $formula = $self -> getFormula();
	my $datatype = $self -> getDatatype();
	my $dbName = $self -> getDBname();	
	if (!$datatype) {
		$datatype = "SMILES";
	}
	if (!$formula) {
		print "Formula is not defined\n";
		return ();
	}
	if (!$database) {
		print "Structure Database path is not defined\n";
		return ();	
	}
	my $jarPath = $self -> getJarPath();
	my $command = "java -jar " . $jarPath . " -QueryStructureDatabase " . $formula . " " . $database . " " . $datatype . " " . $dbName;
	print $command;
    my @value = `$command`;
	return \@value;
}

################################################
## Perform Fragmentation of a SMILE structure ##
################################################
#sub smilesFragmentation {
#	my ($self, $smiles, $mass, $depth) = @_;
#	if (!$smiles) {
#		print "formula is not defined\n";
#		return ();
#	}
#	if (!$mass) {
#		$mass = "0.0";
#	}
#	if (!$depth) {
#		$depth = "2";
#	}
#	my $jarPath = $self -> getJarPath();
#	my $command = "java -jar " . $jarPath . " -FragmentSMILE " . $smiles . " " . $mass . " " . $depth;
#    my @value = `$command`;	
#	return @value;
#}

##################################################
# Calculates the isotope Pattern given a formula #
##################################################
sub AIMPuritySampler {
	my ($self, $formula, $C12File, $C13File, $N15File, $sampledElement, $charge, $ppm) = @_;
	if (!$ppm) {
		$ppm = 5;
	}
	if (!$charge) {
		print "Charge is not defined\n";
		return ();
	}
	if (!$formula) {
		print "Formula is not defined\n";
		return ();	
	}
	if (!$C12File) {
		print "C12 Peak file path is not defined\n";
		return ();	
	}
	if (!$C13File) {
		print "C13 Peak file path is not defined\n";
		return ();	
	}
	if (!$N15File) {
		print "N15 Peak file path is not defined\n";
		return ();	
	}
	if (!$sampledElement) {
		print "The element EX: C12, C13, N15 is not defined\n";
		return ();	
	}
	my $jarPath = $self -> getJarPath();
	my $command = "java -jar " . $jarPath . " -PuritySampler " . $formula . " " . $C12File . " " . $C13File . " " . $N15File . " " . $sampledElement . " " . $charge . " " . $ppm;
    my @value = `$command`;
	return \@value;
}


sub partialCandidate {
	my ($self, $formula, $minMass, $maxMass) = @_;
	my $dtaFile = $self -> getDtaFile();
	my $params = $self -> getParameter();
	my $Ctype = "C13";
	my $minPurity = 0;
	my $maxPurity = 1;
	my $charge = $params -> {'mode'};
	my $ppm = $params -> {'formula_mass_tolerance_searching'};
	if (!$ppm) {
		$ppm = 10;
	}
	if (!$charge) {
		print "Charge is not defined\n";
		return ();
	}
	if (!$formula) {
		print "Formula is not defined\n";
		return ();	
	}
	my $jarPath = $self -> getJarPath();
	my $command = "java -jar /home/tshaw/JumpJarPackage/JumpPackage_20160912.jar" . " -RoughIsotopePurityMeasurementDTA " . $dtaFile . " " . $formula . " " . $Ctype . " " . $charge . " " . $ppm . " " . $minPurity . " " . $maxPurity . " " . $minMass . " " . $maxMass; 
	print $command, "\n";
    my @value = `$command`;
	return \@value;
}

sub getJarPath {
	my ($self) = shift;
	my $Bin = $self -> getBin();
	my $JUMP_PATH = $Bin . "/JumpJarPackage";
	my $jarpath = $JUMP_PATH . "/JumpPackage.jar"; 
    return $jarpath;
}

1;
