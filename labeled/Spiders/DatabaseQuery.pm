#!/usr/local/bin/perl


######### QueryMassFormula ####################################
#                                                             #
#       **************************************************    #  
#       **** QueryMassFormula                         ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2014 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::DatabaseQuery;

use FindBin qw($Bin);
use lib "$Bin";
use POSIX;
no warnings 'uninitialized';


use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.0;

@ISA	 = qw(Exporter);
@EXPORT = qw(setBin getBin setDatabase getDatabase setFormula getFormula setDatatype getDatatype setDBname getDBname setMass getMass set_parameter get_parameter set_dta_file get_dta_file QueryMassWithDatabase  QueryMassWithDatabaseFilter  QueryStructure formula_parser  QueryStructureDatabase  SmileFragmentation  AIMPuritySampler  partial_candidate getJarPath);
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
    };
    bless $self, $class;
    return $self;
}

sub setBin
{
	my ($self,$Bin)=@_;
	$self->{_Bin}=$Bin;
}

sub getBin
{
	my ($self)=@_;
	return $self->{_Bin};	
}

sub setDatabase
{
	my ($self,$database)=@_;
	$self->{_database}=$database;	
}

sub getDatabase
{
	my ($self)=@_;
	return $self->{_database};	
}

sub setFormula
{
	my ($self,$formula)=@_;
	$self->{_formula}=$formula;	
}

sub getFormula
{
	my ($self)=@_;
	return $self->{_formula};	
}

sub setDatatype
{
	my ($self,$datatype)=@_;
	$self->{_datatype}=$datatype;	
}

sub getDatatype
{
	my ($self)=@_;
	return $self->{_datatype};	
}

sub setDBname
{
	my ($self,$DBname)=@_;
	$self->{_dbname}=$DBname;	
}

sub getDBname
{
	my ($self)=@_;
	return $self->{_dbname};	
}

sub setMass
{
	my ($self,$mass)=@_;
	$self->{_mass}=$mass;	
}

sub getMass
{
	my ($self)=@_;
	return $self->{_mass};	
}

sub set_parameter
{
	my ($self,$param)=@_;
	$self->{'_parameter'}=$param;
}


sub get_parameter
{
	my $self=shift;
	return $self->{'_parameter'};
}

sub set_dta_file
{
	my ($self,$dta_file)=@_;
	return $self->{_dta_file}=$dta_file;
}

sub get_dta_file
{
	my $self=shift;
	return $self->{_dta_file};
}


##################################################
# Calculates the isotope Pattern given a formula #
##################################################
sub QueryMassWithDatabase() {
	my ($self,$mass, $ppm, $IsotopeDatabase, $rgdb_check, $dbname) = @_;
	if (!$ppm) {
		$ppm = 1;
	}
	if (!$mass) {
		print "mass is not defined\n";
		return @array=();
	}
	if (!$IsotopeDatabase) {
		print "IsotopeDatabase path is not defined\n";
		return @array=();	
	}
	if (!$rgdb_check) {
		$rgdb_check = "no";
	}
	my $jarpath = $self->getJarPath();
	my $script = "java -jar " . $jarpath . " -queryMassWithPPM " . $mass . " " . $ppm . " " . $IsotopeDatabase . " " . $rgdb_check . " " . $dbname;
    my @value = `$script`;
	return \@value; # @$value
}

##################################################
# Query the formula given a mass                 #
##################################################
sub QueryMassWithDatabaseFilter() {
       my ($self, $mass, $ppm, $IsotopeDatabase, $rgdb_check, $targetDB, $decoyDB) = @_;
       if (!$ppm) {
             $ppm = 1;
       }
       if (!$mass) {
             print "mass is not defined\n";
             return @array=();
       }
       if (!$IsotopeDatabase) {
             print "IsotopeDatabase path is not defined\n";
             return @array=();   
       }
       if (!$rgdb_check) {
             $rgdb_check = "no";
       }
       if (!$targetDB) {
             print "TargetDB term is not defined\n";
             return @array=();
       }
       if (!$decoyDB) {
             print "DecoyDB term is not defined\n";
             return @array=();
       }
       my $jarpath = $self->getJarPath();
       my $script = "java -jar " . $jarpath . " -QueryMassWithPPMFilter " . $mass . " " . $ppm . " " . $IsotopeDatabase . " " . $rgdb_check . " " . $targetDB . " " . $decoyDB;
    my @value = `$script`;
    
    my %all;
    my %target;
    my %decoy;
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
    my @all_array;
    for (keys %all) {
        push @all_array, $_;
    }
    my @target_array;
    for (keys %target) {
        push @target_array, $_;
    }
    my @decoy_array;
    for (keys %decoy) {
        push @decoy_array, $_;
    }
    
    return (\@all_array, \@target_array, \@decoy_array);
}


sub QueryStructure
{
	my ($self) = @_;
	my @str = ();
	my $database = $self->getDatabase();
	my $mass = $self->getMass();
	my $formula = $self->getFormula();
	my $formula_file = $formula;
	$formula_file =~ s/\d+//g;

	my ($mass_int,$mass_dec) = split(/\./,$mass);
	my $structure_file = "$database/$mass_int/$formula_file" . ".txt";
	
	if(-e($structure_file))
	{
		open(FILE,$structure_file);
		while(<FILE>)
		{
			my @data = split(/\t/,$_);
			next if(!defined($formula));
			next if(!defined($data[1]));
			
			my $element1 = $self->formula_parser($formula);
			my $element2 = $self->formula_parser($data[1]);			

			my $formula1="";
			my $formula2="";
			foreach my $key (sort keys %$element1)
			{
				my $elm = $key . $element1->{$key};
				$formula1 .= $elm;
			}


			foreach my $key (sort keys %$element2)
			{
				my $elm = $key . $element2->{$key};
				$formula2 .= $elm;
			}
			
			if($formula1 eq $formula2)
			{			
				my $smile_iupac = join("\t",$data[2],@data[4..5]);
				push(@str,$smile_iupac);
			}
		}
	}
	return \@str;
}

sub formula_parser {
    my ($self,$formula) = @_;
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


##################################################
# Query the structure given a formula #
# $datatype = INCHIKEY;INCHISTR;SMILE,IUPAC,GENERALNAME, (multiple selection)
# $dbname = PUBCHEM,HMDB,KEGG
##################################################
sub QueryStructureDatabase() {
#	my ($self,$formula, $database, $datatype) = @_;
	$self=shift;
	my $database = $self->getDatabase();
	my $formula = $self->getFormula();
	my $datatype = $self->getDatatype();
	my $db_name = $self->getDBname();
	
	
	if (!$datatype) {
		$datatype = "SMILES";
	}
	if (!$formula) {
		print "formula is not defined\n";
		return @array=();
	}
	if (!$database) {
		print "Structure Database path is not defined\n";
		return @array=();	
	}
	my $jarpath = $self->getJarPath();
	my $script = "java -jar " . $jarpath . " -QueryStructureDatabase " . $formula . " " . $database . " " . $datatype . " " . $db_name;
	print $script;
    my @value = `$script`;
	return \@value;
}
##################################################
# Perform Fragmentation of a SMILE structure #
##################################################
sub SmileFragmentation() {
	my ($self,$smile, $mass, $depth) = @_;
	if (!$smile) {
		print "formula is not defined\n";
		return @array=();
	}
	if (!$mass) {
		$mass = "0.0";
	}
	if (!$depth) {
		$depth = "2";
	}
	my $jarpath = $self->getJarPath();
	my $script = "java -jar " . $jarpath . " -FragmentSMILE " . $smile . " " . $mass . " " . $depth;
    my @value = `$script`;
	
	return @value;
}
##################################################
# Calculates the isotope Pattern given a formula #
##################################################
sub AIMPuritySampler() {
	my ($self,$formula, $C12_File, $C13_File, $N15_File, $sampled_element, $charge, $ppm) = @_;
	if (!$ppm) {
		$ppm = 5;
	}
	if (!$charge) {
		print "charge is not defined\n";
		return @array=();
	}
	if (!$formula) {
		print "formula is not defined\n";
		return @array=();	
	}
	if (!$C12_File) {
		print "C12 Peak file path is not defined\n";
		return @array=();	
	}
	if (!$C13_File) {
		print "C13 Peak file path is not defined\n";
		return @array=();	
	}
	if (!$N15_File) {
		print "N15 Peak file path is not defined\n";
		return @array=();	
	}
	if (!$sampled_element) {
		print "The element EX: C12, C13, N15 is not defined\n";
		return @array=();	
	}
	my $jarpath = $self->getJarPath();
	my $script = "java -jar " . $jarpath . " -PuritySampler " . $formula . " " . $C12_File . " " . $C13_File . " " . $N15_File . " " . $sampled_element . " " . $charge . " " . $ppm;
    my @value = `$script`;

	return \@value;
}


sub partial_candidate 
{
	my ($self,$formula, $min_mass, $max_mass) = @_;
	my $dta_file = $self->get_dta_file();
	my $params = $self->get_parameter();
	
	my $C_type = "C13";
	my $min_purity = 0;
	my $max_purity = 1;
	my $charge = $params->{'mode'};
	my $ppm = $params->{'formula_mass_tolerance_searching'};
	
	if (!$ppm) {
		$ppm = 10;
	}
	if (!$charge) {
		print "charge is not defined\n";
		return @array=();
	}
	if (!$formula) {
		print "formula is not defined\n";
		return @array=();	
	}

	my $jarpath = $self->getJarPath();

	
	my $script = "java -jar /home/tshaw/JumpJarPackage/JumpPackage_20160912.jar" . " -RoughIsotopePurityMeasurementDTA " . $dta_file . " " . $formula . " " . $C_type . " " . $charge . " " . $ppm . " " . $min_purity . " " . $max_purity . " " . $min_mass . " " . $max_mass; 
	print $script,"\n";
    my @value = `$script`;

	return \@value;
}


sub getJarPath
{
	my ($self) = shift;
	my $Bin = $self->getBin();
	my $JUMP_PATH = $Bin . "/JumpJarPackage";
	my $jarpath = $JUMP_PATH . "/JumpPackage.jar"; 
#	my $jarpath = $JUMP_PATH . "/JumpPackage_20160901.jar"; 

    return $jarpath;
}

1;
