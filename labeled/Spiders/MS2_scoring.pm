#!/usr/bin/perl 

######### Pairing #############################################
#                                                             #
#       **************************************************    #  
#       **** Pairing                    	          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2014 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::MS2_scoring;
use Spiders::Hypergeometric;
use Spiders::WilcoxonRankSum;

use warnings;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.0;

@ISA	 = qw(Exporter);
@EXPORT = qw(compare_theoritical_experiment wilcox_test get_pvalue set_exp_mz get_exp_mz get_exp_mz_hash set_parameter get_parameter);
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
    };
    bless $self, $class;
    return $self;
}

sub compare_theoritical_experiment
{
	my ($self,$exp_mz,$theor_mz) = @_;
	my $parameter = $self->get_parameter();
	my $matched = 0;
	my $tolerance_input = $parameter->{'frag_mass_tolerance'};
	my $tolerance_unit = 1;
	if(defined($parameter->{'frag_mass_tolerance_unit'}))
	{                        
		$tolerance_unit = $parameter->{'frag_mass_tolerance_unit'};
	}

	
#### 6/28/2013 Fixed a bug: count multiple times because of similar theoretical mass generated from 1 charge and 2 charge	
	my %matched_hash;

	$tolerance = $tolerance_input;
### version 11.0.2 to boot the speed

	for (my $i=0;$i<=$#{$theor_mz};$i++)
	{	
	
		my $theor_mz_value = $theor_mz->[$i];
		my $theor_mz_value_int = int($theor_mz_value);

		if($tolerance_unit == 2)
		{
			$tolerance = $tolerance_input * $theor_mz_value_int / 1000000;
		}
		my $theor_mz_value_min = int($theor_mz_value-$tolerance);
		my $theor_mz_value_max = int($theor_mz_value+$tolerance)+1;
	
		for (my $theor_mz_value_i = $theor_mz_value_min; $theor_mz_value_i < $theor_mz_value_max; $theor_mz_value_i++)
		{
			if(defined($exp_mz->{$theor_mz_value_i}))
			{
				foreach my $exp_mz_value (keys %{$exp_mz->{$theor_mz_value_i}})
				{
				
					if(($exp_mz_value+$tolerance)>$theor_mz_value and ($exp_mz_value-$tolerance)<$theor_mz_value)
					{
						next if($matched_hash{$exp_mz_value});

						$matched++;
#						print "matched:$theor_mz_value\t$exp_mz_value\n";
						$matched_hash{$exp_mz_value}=1;
					}

				}
			}
		}			
	}

	$matched = scalar keys (%matched_hash);
#	print "matched $matched peaks\n";
	return ($matched,\%matched_hash);
}

sub wilcox_test
{
	my ($self,$matched_array_ref,$total_array_ref)=@_;
    use Spiders::WilcoxonRankSum;

    my $wilcox_test = Spiders::WilcoxonRankSum->new();

	my @matched_array=@$matched_array_ref;
	my @total_array=@$total_array_ref;
	my %matched_hash=map{$_=>1} @matched_array;
	my @unmatched_array=grep(!defined $matched_hash{$_}, @total_array);
	my $prob = 1;
	if(scalar(@matched_array)>0)
	{
		$wilcox_test->load_data(\@matched_array, \@unmatched_array);
		$prob = $wilcox_test->probability();
	}
    my $pf = sprintf '%f', $prob; # prints 0.091022
	return $pf;	
}

sub get_pvalue
{
	my ($self,$matched,$total)=@_;
######### use different method for matching score #############
	my $parameter = $self->get_parameter();
	
###### Version 1.12 fix a bug: using frag tolerance instead of mass tolerance	
	my $mass_tolerance = $parameter->{'frag_mass_tolerance'};
	
	my $tolerance_unit = 1;
	if(defined($parameter->{'frag_mass_tolerance_unit'}))
	{
		$tolerance_unit = $parameter->{'frag_mass_tolerance_unit'};
	}
	
	my $exp_mz = $self->get_exp_mz();
	
	if($tolerance_unit == 2)
	{
### use the average mass to convert the unit	
		$mass_tolerance = $mass_tolerance *($exp_mz->[$#$exp_mz] + $exp_mz->[0]) / 1000000;
	}			
	
	my $total_number = int($exp_mz->[$#$exp_mz] - $exp_mz->[0])/($mass_tolerance*2);


	my $exp_mass_num = scalar (@$exp_mz);
	
	my $hyper = new Spiders::Hypergeometric();

	my $log_peptide_pvalue = 1;

	my $peptide_pvalue=$hyper->Hypergeometric($total_number,$exp_mass_num,$total,$matched);	
	
	return ($peptide_pvalue);
}

sub set_exp_mz
{
	my ($self,$exp_mz)=@_;
	$self->{'_exp_mz'} = $exp_mz;
}

sub get_exp_mz
{
	my $self=shift;
	return $self->{'_exp_mz'};
}

sub get_exp_mz_hash
{
	my $self=shift;
	my $exp_mz = $self->{'_exp_mz'};
	my %exp_hash=();
	foreach (@$exp_mz)
	{
		my $int_mz = int($_);
		$exp_hash{$int_mz}{$_}=1;
	}
	return \%exp_hash;
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

1;







