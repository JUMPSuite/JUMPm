#!/usr/bin/perl

############################################################# 
#                                                           #
#       **************************************************  #    
#       **** Decharge program for MS1 and MS2         ****  #    
#       ****                                          ****  #    
#       ****Copyright (C) 2012 - Xusheng Wang        ****  #    
#       ****all rights reserved.                      ****  #    
#       ****xusheng.wang@stjude.org                   ****  #    
#       ****                                          ****  #    
#       ****                                          ****  #    
#       **************************************************  #    
#############################################################

package Spiders::Decharge;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT = qw(set_C_value get_C_value set_H_value get_H_value set_parameter get_parameter set_dta_path get_dta_path find_charge_mono decharge MS1_deisotope get_isotopic_distribution get_intensity_ratio find_charge); 

sub new{
	my ($class,%arg)=@_;
    my $self = {
		_C_value =>undef,
		_H_value=>undef,
    };
    bless $self, $class;
	return $self;
}

sub set_C_value
{
	my ($self,$c_value)=@_;
	$self->{_C_value}=$c_value;	
}

sub get_C_value
{
	my $self=shift;
	if(!defined($self->{_C_value}))
	{
		$self->{_C_value}=1.00335;
	}
	return $self->{_C_value};
}

sub set_H_value
{
	my ($self,$c_value)=@_;
	$self->{_H_value}=$c_value;	
}

sub get_H_value
{
	my $self=shift;
	if(!defined($self->{_H_value}))
	{
		$self->{_H_value}=1.007276466812;
	}
	return $self->{_H_value};
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

sub set_dta_path
{
	my ($self,$dir)=@_;
	$self->{'_dta_path'} = $dir;
}

sub get_dta_path
{
	my $self=shift;
	return $self->{'_dta_path'};
}

sub find_charge_mono
{

}


sub decharge{
	my ($self,  $mshash, $msmshash, $mzarray) = @_;
	my $H = $self->get_H_value();
	
	my $parameter = $self->get_parameter();
### keep consistent with parameter file	
	my $iso_window = $parameter->{'isolation_window'};
	# Decharge

	
	
	
	foreach my $specscan  (keys %$msmshash)
	{

		my $sourcemz = 	$msmshash->{$specscan}->{'prec_mz'};
		
		my $survey = $msmshash->{$specscan}->{'survey'};
		next if (!defined($sourcemz));

  # Use original scan to begin deisotoping
		my $charge = $self->find_charge($mshash, $mzarray, $survey, $sourcemz);

		my ($low, $high) = ($sourcemz-$iso_window/2, $sourcemz+$iso_window/2);
	  
		my @lowarray = split('\.', $low);
		my @higharray = split('\.', $high);
		my ($lowint, $highint) = ($lowarray[0], $higharray[0]);

	  # Create a hash with all mz values between lowint and highint
		my %mzhash;

	   #    print "$strongest_mz low = $lowint high = $highint\n";
		for (my $i = $lowint; $i<=$highint; $i++)
		{
			next if (!defined(%{$$mzarray[$survey][$i]}));
			while (my ($key, $value) = each %{$$mzarray[$survey][$i]})
			{

				if($key>=$low and $key<=$high)
				{
					$mzhash{$key} = $value;	
				}		
			}
		}
		my $strongest_mz = $sourcemz;
		$mzhash{$strongest_mz} = 0;		
		foreach my $mz (keys %mzhash)
		{
			if($mzhash{$mz}>$mzhash{$strongest_mz})
			{
				$strongest_mz = $mz;
			}
		}		
		
		my $updated_mz = $strongest_mz; 	
		$msmshash->{$specscan}->{'prec_MH'} = (sprintf("%.6f", $updated_mz) - $H) * $charge + $H;		
	}	

}

sub MS1_deisotope
{
	my ($self,$select_mz,$charge,$mz_hash)=@_;

	my $C = $self->get_C_value();

	my $parameter = $self->get_parameter();
### keep consistent with parameter file	
	my $tolerance = $parameter->{'deisotope_ppm'};

	my $prec_mz = $select_mz;
	my $flag = 1;
	my $max_loop = $self->get_isotopic_distribution($select_mz*$charge);
	
	for(my $search_loop=0; $search_loop < $max_loop; $search_loop++)
	{
		my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz - ($C/$charge)*$search_loop - $tolerance, $select_mz - ($C/$charge)*$search_loop + $tolerance); 
		foreach my $mz ( keys %{$mz_hash})
		{
	# search the previous peak ( only one peak) 

			if($mz>$previous_peak_mz_low && $mz<$previous_peak_mz_high)
			{

				if(defined ($mz_hash->{$mz}) && defined ($mz_hash->{$select_mz}) )
				{
					my $intensity_ratio = $self->get_intensity_ratio($select_mz*$charge,$search_loop);

					if(($mz_hash->{$mz} / $mz_hash->{$select_mz}) > $intensity_ratio)
					{
						$prec_mz = $mz;
					}
				}
			}
		}
	}	
	return $prec_mz;
}

sub get_isotopic_distribution
{
	my ($self, $mz) = @_;
	my $loop = 0;
	if($mz<1500)
	{
		$loop = 1;
	}
	elsif($mz<3000 && $loop<=2)
	{
		$loop = 2;
	}
	elsif($mz<4500 && $loop<=4)
	{
		$loop = 4;
	}
	elsif($mz<6000 && $loop<=6)
	{
		$loop = 6;
	}	
	return $loop;	
	
}

sub get_intensity_ratio
{
	my ($self,$mz,$loop) = @_;
	my $ratio=2;
### if the mass <1500, there is no isotopic peaks preceeding the monoisotope	
	if($mz<1500)
	{
		$ratio = 0.8;
	}
	elsif($mz<3000 && $loop<=2)
	{
		if($loop == 1)
		{
			$ratio = 0.4;
		}
		elsif($loop == 2)
		{
			$ratio = 0.3;
		}
	}
	elsif($mz<4500 && $loop<=4)
	{
		if($loop == 1)
		{
			$ratio = 0.6;
		}
		elsif($loop == 2)
		{
			$ratio = 0.2;
		}
		elsif($loop == 3)
		{
			$ratio = 0.1;
		}
		elsif($loop == 4)
		{
			$ratio = 0.1;
		}
	}
	elsif($mz<6000 && $loop<=6)
	{
		if($loop == 1)
		{
			$ratio = 0.6;
		}
		elsif($loop == 2)
		{
			$ratio = 0.3;
		}
		elsif($loop == 3)
		{
			$ratio = 0.1;
		}
		elsif($loop == 4)
		{
			$ratio = 0.1;
		}
		elsif($loop == 5)
		{
			$ratio = 0.1;
		}
		elsif($loop == 6)
		{
			$ratio = 0.1;
		}
	}	
	return $ratio;
}



sub find_charge{
	my ($self, $mshash, $mzarray, $scan, $origmz) = @_;

	my $charge=1;

	my $parameter = $self->get_parameter();
	my $intrappm = $parameter->{'decharge_ppm'};

	my $C = $self->get_C_value();

	my $maxcharge = 5;

	my $strongest_mz = $origmz;
	my $strongest_scan = $scan;
####### Version 1.15: define mass region: First step to find whether there is any charge between 2 to 5 using $C/1.8
#################### 
#############
	my ($low, $high) = ($strongest_mz-$C-($intrappm/1000000)*$strongest_mz, $strongest_mz+$C+($intrappm/1000000)*$strongest_mz);
  
	my @lowarray = split('\.', $low);
	my @higharray = split('\.', $high);
	my ($lowint, $highint) = ($lowarray[0], $higharray[0]);

  # Create a hash with all mz values between lowint and highint
	my %mzhash;

	for (my $i = $lowint; $i<=$highint; $i++){
		next if (!defined(%{$$mzarray[$strongest_scan][$i]}));
		while (my ($key, $value) = each %{$$mzarray[$strongest_scan][$i]}){
###### selected the peaks within the exact window, not the expanded window	
			if($key>=$low and $key<=$high)
			{
				$mzhash{$key} = $value;	
			}
	
		}
	}


	foreach my $mz (keys %mzhash)
	{
		next if(defined($mzhash{$strongest_mz}));
		next if(defined($mzhash{$mz}));		
		if($mzhash{$mz}>$mzhash{$strongest_mz})
		{
			$strongest_mz = $mz;
		}
	}
######### get the strongest, then selected again
	
	($low, $high) = ($strongest_mz-$C-($intrappm/1000000)*$strongest_mz, $strongest_mz+$C+($intrappm/1000000)*$strongest_mz);
	@lowarray = split('\.', $low);
	@higharray = split('\.', $high);	
	($lowint, $highint) = ($lowarray[0], $higharray[0]);

	undef %mzhash;
	for (my $i = $lowint; $i<=$highint; $i++){
		next if (!defined(%{$$mzarray[$strongest_scan][$i]}));
		while (my ($key, $value) = each %{$$mzarray[$strongest_scan][$i]})
		{
				$mzhash{$key} = $value;	
		}
	}	
	
  # Create a hash of only the peaks between low and high
	my %found;
	while (my ($mz, $intensity) = each  %mzhash){
		next if($mz<$low || $mz>$high);
		$found{$strongest_scan}{$mz} = $intensity;
	}

	my $diffsum=0;

  # Get mz with the highest intensity
	for my $mz (sort {$found{$strongest_scan}{$b}<=>$found{$strongest_scan}{$a}} keys %{$found{$strongest_scan}}){
		next if ($mz==$strongest_mz);
		my $diff = 1/abs($mz-$strongest_mz);
		my $round_diff = sprintf("%.0f", $diff);
		next if ($round_diff==0);
		next if ($round_diff > $maxcharge);
		my $var = abs(abs($mz-$strongest_mz)-($C/$round_diff));
		next if ($var > ($intrappm/1000000)*$strongest_mz);
                
		$diffsum += ($var*1000000)/$mz;
		$charge = $round_diff;

		last;
	}
	
	
	return $charge;
	
}



1;
