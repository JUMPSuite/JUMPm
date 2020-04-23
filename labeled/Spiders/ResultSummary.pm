#!/usr/bin/perl

######### Result summary ######################################
#                                                             #
#       **************************************************    #  
#       **** Search Structure for Metabolite          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::ResultSummary;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.03;

@ISA	 = qw(Exporter);
@EXPORT      = ();
 
set_smout_folder($ARGV[0]); 
Readsmout(); 

sub new{
    my ($class,%arg) = @_;
    my $self = {
    };
    bless $self, $class;
    return $self;
}

sub set_smout_folder
{
	my ($self,$folder) = @_;
	$self->{'smout_folder'} = $folder;	
}

sub get_smout_folder
{
	my $self = shift;
	return $self->{'smout_folder'};
}

sub Readsmout
{
	my $self = @_;
	my $smout_folder = $self->get_smout_folder();
	my @smout_files = glob("$smout_folder/*.smout");
	foreach my $file (@smout_files)
	{
		open(SMOUT,$file) || die "can not open the file: $file";
		while(<SMOUT>)
		{
			$last_line = $_;
		}
		my @data = split(/\t/,$last_line);
		print $data[0],"\n";
	}
}