#!/usr/local/bin/perl


######### QueryMassFormula #################################################
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


use POSIX;
no warnings 'uninitialized';


use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.0;

@ISA	 = qw(Exporter);
@EXPORT      = ();
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
    };
    bless $self, $class;
    return $self;
}

##################################################
# Calculates the isotope Pattern given a formula #
##################################################
sub QueryMassWithDatabase() {
	my ($mass, $ppm, $IsotopeDatabase, $rgdb_check) = @_;
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
	my $jarpath = getJarPath();
	my $script = "java -jar " . $jarpath . " -queryMassWithPPM " . $mass . " " . $ppm . " " . $IsotopeDatabase . " " . $rgdb_check;
    my @value = `$script`;
	return \@value; # @$value
}
##################################################
# Query the structure given a formula #
##################################################
sub QueryStructureDatabase() {
	my ($formula, $database, $datatype) = @_;
	if (!$datatype) {
		$datatype = "SMILE";
	}
	if (!$formula) {
		print "formula is not defined\n";
		return @array=();
	}
	if (!$database) {
		print "Structure Database path is not defined\n";
		return @array=();	
	}
	my $jarpath = getJarPath();
	my $script = "java -jar " . $jarpath . " -QueryStructureDatabase " . $formula . " " . $database . " " . $datatype;
    my @value = `$script`;
	return @value;
}
##################################################
# Perform Fragmentation of a SMILE structure #
##################################################
sub SmileFragmentation() {
	my ($smile, $mass, $depth) = @_;
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
	my $jarpath = getJarPath();
	my $script = "java -jar " . $jarpath . " -FragmentSMILE " . $smile . " " . $mass . " " . $depth;
    my @value = `$script`;
	return @value;
}
##################################################
# Calculates the isotope Pattern given a formula #
##################################################
sub AIMPuritySampler() {
	my ($formula, $C12_File, $C13_File, $N15_File, $sampled_element, $charge, $ppm) = @_;
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
	my $jarpath = getJarPath();
	my $script = "java -jar " . $jarpath . " -PuritySampler " . $formula . " " . $C12_File . " " . $C13_File . " " . $N15_File . " " . $sampled_element . " " . $charge . " " . $ppm;
    my @value = `$script`;
	return @value;
}


sub getJarPath{
	my $JUMP_PATH = "/spiders/home/xwang4/scripts/JUMPm/JumpJarPackage";
	print "JUMPJARPATH is: $JUMP_PATH\n";
	my $jarpath = $JUMP_PATH . "/JumpPackage_0.5.jar"; # this might only work for version 0.1 if not changed
	return $jarpath;	
}

1;
