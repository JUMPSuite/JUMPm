#!/usr/bin/perl


######### Deisotope ##########################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2012 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################
# this value is used to remove the precursor ion

package Spiders::Deisotope;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);
use File::Basename;

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT = qw(set_C_value get_C_value set_H_value get_H_value set_dta get_dta set_parameter get_parameter MS2_deisotope MS1_deisotope calculate_signal_noise_ratio median remove_prec sort_hash define_charge deisotope get_isotopic_peaks_mass_error print_mass_error changeMH deisotope_charge_1);  

sub new{
	my ($class,%arg)=@_;
    my $self = {
		_C_value =>undef,
		_H_value=>undef,
		_mass_error=>undef,
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


sub set_dta
{
	my ($self,$dta)=@_;
	$self->{'_dta_file'} = $dta;
	return $self->{'_dta_file'};
}


sub get_dta
{
	my $self=shift;
	return $self->{'_dta_file'};
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

sub MS2_deisotope
{
	my $self=shift;
	my $H=$self->get_H_value();

	my $dta = $self->get_dta();
	my $parameter=$self->get_parameter();

# define the mz hash for storing the mz and intensity
    my %mz_hash;
# open each dta file to load the data 

    open(DTAFILE,$dta) || die "can not open the dta file: $dta";
# get the precusor mz and charge 
    my $prec_mz_charge = <DTAFILE>;
    my ($prec_mass,$prec_charge) = split(/\s+/,$prec_mz_charge);

    while(<DTAFILE>)
    {
# get the mz and intensity and save in the mz_hash
		my @data =split(/\s+/,$_);
		$mz_hash{$data[0]} = $data[1];
	}
# sort the mz according to the intensity

    my @mz_array;
	my %prec_mz_hash;
# remove un-fragmented ion mass (3-unit), save the precursor value in %prec_mz_hash
	my $prec_mz = (($prec_mass-$H)/$prec_charge)+$H;
	$self->remove_prec(\%mz_hash,$prec_mz,\%prec_mz_hash);

## sort mz according to the intensity


#	$self->sort_hash(\%mz_hash,\@mz_array);

	my %charge_hash;
#decharge and deiotope for each peaks 
#	foreach my $mz (@mz_array)
	foreach my $mz (reverse sort {$mz_hash{$a}<=>$mz_hash{$b}} keys %mz_hash)
	{

		next if (!defined($mz_hash{$mz}));
# find charge and isotopes
		my $charge = $self->define_charge(\%mz_hash,$parameter->{'ppm'},$mz,$prec_charge);
## define the isotopic peaks using the charge, tolerance and orignal peaks
		next if ($charge == 0);
		if(defined($charge))
		{
			$self->deisotope(\%mz_hash,$charge,$mz,$parameter->{'ppm'});
			$charge_hash{$mz}=$charge;
		}
	}

### change the mass if the charge is larger than 1
#	print "debuging........................\n";
	foreach my $mz (keys %mz_hash)
	{
#		next if (!defined $charge_hash{$mz});
		if(defined $charge_hash{$mz} && $charge_hash{$mz} >1)
		{
				$self->changeMH(\%mz_hash,$charge_hash{$mz},$mz);
		}

	}

	my $ms2_signal_noise_ratio = $self->calculate_signal_noise_ratio($prec_mass,$prec_charge,\%mz_hash);
	
	open(OUT,">$dta");
	print OUT $prec_mz_charge;

	foreach my $mz (sort {$a<=>$b} keys %mz_hash)
	{
		if(defined($mz_hash{$mz}))
		{
			print OUT $mz," ",$mz_hash{$mz},"\n";
		}
	}
	close(OUT);
	return $ms2_signal_noise_ratio;
}

sub MS1_deisotope
{
	my $self=shift;
	my $H=$self->get_H_value();

	my $dta = $self->get_dta();
	my $parameter=$self->get_parameter();

# define the mz hash for storing the mz and intensity
    my %mz_hash;
	my %int_ratio_hash;
# open each dta file to load the data 

    open(DTAFILE,$dta) || die "can not open the dta file: $dta";
# get the precusor mz and charge 
#    my $prec_mz_charge = <DTAFILE>;


    while(<DTAFILE>)
    {
# get the mz and intensity and save in the mz_hash
		my @data =split(/\s+/,$_);
		$mz_hash{$data[0]} = $data[1];

	}
# sort the mz according to the intensity

    my @mz_array;
	my %prec_mz_hash;
	my %peak_C;


## sort mz according to the intensity


#	$self->sort_hash(\%mz_hash,\@mz_array);

	my %charge_hash;
#decharge and deiotope for each peaks 
#	foreach my $mz (@mz_array)

	foreach my $mz (reverse sort {$mz_hash{$a}<=>$mz_hash{$b}} keys %mz_hash)
	{
		next if (!defined($mz_hash{$mz}));
		next if(defined($charge_hash{$mz}));
# find charge and isotopes
		my $charge=0;
		my $int_ratio=0;
		my $other_peak_mz = 0;
		my $peak_C_flag="";
# $other_peak_mz for the right/left side peak
		($charge,$int_ratio,$other_peak_mz,$peak_C_flag) = $self->define_charge(\%mz_hash,$parameter->{'decharge_ppm'},$mz,6);
		$int_ratio_hash{$mz} = $int_ratio;
## define the isotopic peaks using the charge, tolerance and orignal peaks
		next if ($charge == 0);
		if(defined($charge))
		{
			$charge_hash{$mz}=$charge;
			$charge_hash{$other_peak_mz}=$charge;			
			$peak_C{$mz}=$peak_C_flag;			
		}		
	}
	
	foreach my $mz (sort {$mz_hash{$b}<=>$mz_hash{$a}} keys %mz_hash)
	{
		if(defined $charge_hash{$mz} && $charge_hash{$mz} >1)
		{				
				my $new_mz = ($mz-$H)*$charge_hash{$mz}+$H;
				$peak_C{$new_mz}=$peak_C{$mz};				
				$self->changeMH(\%mz_hash,$charge_hash{$mz},$mz);
		}
	}
	
	return (\%mz_hash,\%int_ratio_hash,\%peak_C,\%charge_hash);
}


sub calculate_signal_noise_ratio
{
	my ($self,$prec_mass,$charge,$mz_hash) = \@_;
	my $charge_value = ($charge >= 2) ? 2 : 1; 
	my $number_ions_select = int(($prec_mass / 118)*0.4) * $charge_value;
	
	my @signal_intensity_array=();
	my @noise_intensity_array=();
	my $i=0;
	if(!(defined(%$mz_hash)))
	{
		return "N/A";
	}
	elsif(scalar(keys %$mz_hash)<2)
	{
		return 0;
	}
	foreach my $mz (reverse sort {$mz_hash->{$a}<=>$mz_hash->{$b}} keys %$mz_hash)
	{
		next if ($mz<150);
		if($i<$number_ions_select)
		{
			push(@signal_intensity_array,$mz_hash->{$mz});
		}
		else
		{
			push(@noise_intensity_array,$mz_hash->{$mz});			
		}
		$i++;
	}
	if(scalar (@signal_intensity_array) < 2 || scalar (@noise_intensity_array)<2)
	{
		return "N/A";
	}
	
	my $signal_median = $self->median(\@signal_intensity_array);
	my $noise_median = $self->median(\@noise_intensity_array);
	if($noise_median ==0)
	{
		return "N/A";
	}
	my $signal_noise_ratio = $signal_median / $noise_median;
	return $signal_noise_ratio;
}

sub median
{
	my ($self,$array)=@_;
    my @vals = sort {$a <=> $b} @$array;
    my $len = @vals;
	return 0 if($len == 0);
	
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
		if(!defined($vals[int($len/2)]))
		{
			return $vals[int($len/2)-1];
		}
		elsif(!defined($vals[int($len/2)-1]))
		{
			return 0;
		}
		else
		{
			return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
		}
    }
}

######### Remove un-fragmented ions #################################
# remove the precursor mass according the window defined by user
# normally, the window size is 3
# save the precursor mass in %prec_mz_hash
# remove the keys and values from %mz_hash

sub remove_prec
{
	my ($self, $mz_hash,$prec_mz,$prec_mz_hash) = @_;
	my $C = $self->get_C_value();
	my $parameter=$self->get_parameter();

	my ($low, $high) = ($prec_mz- $parameter->{'prec_window'}/2*$C, $prec_mz + $parameter->{'prec_window'}/2*$C);
	while (my ($mz, $intensity) = each  %$mz_hash)
	{
		if($mz>$low && $mz<$high)
		{
			$$prec_mz_hash{$mz} = $intensity;
			delete $$mz_hash{$mz};
		}
	}
}

sub sort_hash
{
	my $self=shift;
	my ($mz_hash,$mz_array) = @_;

	my @mz_array=();
	foreach  my $mz (reverse sort {$$mz_hash{$a}<=>$$mz_hash{$b}} keys %$mz_hash)
	{
		push (@$mz_array,$mz);
	}
}


# M=1/M 

########## within one unit ################################################
# the largest charge can be obtained from precursor ion charge
# if the charge between two peaks > precusor ion charge, then go to next peak
#############################################################################



sub define_charge
{
### input #######
# $mz_hash is the reference of %mz_hash
# $mz is the specific mass that is used to define the window for it
# $ppm is the mass tolerence used for defining window
#################
	my $self=shift;
	my ($mz_hash,$ppm,$orign_mz,$prec_charge) = @_;
	my $parameter=$self->get_parameter();
	my $C = $self->get_C_value();

## The max charge can not larger than that of precursor ion
	my $maxcharge = $prec_charge;
	my $charge;
	my $int_ratio=0;
######## only for matebolite, change it to narrow window to exclude one charge and define charge for all peaks
	my ($low, $high) = ($orign_mz-$C-($parameter->{'decharge_ppm'}/1000000)*$orign_mz, $orign_mz+$C+($parameter->{'decharge_ppm'}/1000000)*$orign_mz);
	my @lowarray = split('\.', $low);
	my @higharray = split('\.', $high);
	my ($lowint, $highint) = ($lowarray[0], $higharray[0]);

	my %found;
	my %charge_found;
	while (my ($mz, $intensity) = each  %$mz_hash)
	{
		next if($mz<$low || $mz>$high); 
	    $found{$mz} = $intensity;
	}
  # Get mz with the highest intensity
	my $strongest_mz = $orign_mz;
	my $other_peak_mz = 0;
	my $peak_C="";
	for my $mz (sort {$found{$b}<=>$found{$a}} keys %found)
	{
		next if ($found{$mz}>=$found{$strongest_mz});

		my $diff = 1/abs($mz-$strongest_mz);
		my $round_diff = sprintf("%.0f", $diff);
		next if ($round_diff > $maxcharge);
		my $var = abs(abs($mz-$strongest_mz)-($C/$round_diff));

		next if ($var > $parameter->{'decharge_ppm'}/1000000*$strongest_mz);
		$int_ratio = ($found{$mz}/$found{$strongest_mz})*100;
		my $diffsum += ($var*1000000)/$mz;
		$charge = $round_diff;
		$other_peak_mz = $mz;
		if($mz<$strongest_mz)
		{
			$peak_C = "C13";
		}
		else
		{
			$peak_C = "C12";			
		}
		last;
	}
	if(!defined($charge))
	{
		$charge=0;		
	}

	
	
	return ($charge,$int_ratio,$other_peak_mz,$peak_C);

}

sub deisotope
{
	my $self=shift;
	my $parameter=$self->get_parameter();
	my $C = $self->get_C_value();
	my $delete_mz;
	my ($mz_hash,$charge,$select_mz) = @_;


# search the following peaks 

	my $search_loop = 1;
	my $flag=1;
	my $previous_int = $mz_hash->{$select_mz};
	my $selected_int = $mz_hash->{$select_mz};
	while($search_loop && $flag)
	{
		
		my ($peak_mz_low,$peak_mz_high) =  ($select_mz + ($C/$charge)*$search_loop-($parameter->{'deisotope_ppm'}/1000000)*$select_mz, $select_mz + ($C/$charge)*$search_loop +($parameter->{'deisotope_ppm'}/1000000)*$select_mz);

		
		$flag=0;
		foreach my $mz (keys %$mz_hash)
		{

			next if($mz<100);

			my $previous_mz;
			if($mz>$peak_mz_low && $mz<$peak_mz_high)
			{
				$flag=1;
				
				if(($mz_hash->{$mz} / $previous_int) > 0.3)
				{
					$previous_int = $mz_hash->{$mz};
					$previous_mz = $mz;
					

					if(!defined($mz_hash->{$select_mz}))
					{
						$search_loop++;
						next;
					}
	###### always keep the strongest peaks....... (only for metabolite detection)				
					if($mz_hash->{$mz}>$mz_hash->{$select_mz})
					{					
						$mz_hash->{$mz} += $mz_hash->{$select_mz};		
						$search_loop++;

						$delete_mz = $select_mz;
					
					}
					else
					{
						$mz_hash->{$select_mz} += $mz_hash->{$mz};
						$search_loop++;
						$delete_mz = $mz;						

					}
				}				
				else
				{
					$search_loop=0;
					$previous_mz=0;
				}
			}
		
		}
		
	}
	return $delete_mz;
		
}       

sub get_isotopic_peaks_mass_error
{
	my $self = shift;
	return $self->{'_mass_error'};
}

sub print_mass_error
{
	my ($self,$filename) = @_;
	my $mass_error = $self->get_isotopic_peaks_mass_error();
	open(FILE,">$filename") || die "can not open the file";
	foreach my $mz (keys %$mass_error)
	{
		
		print FILE $mz,"\t",$mass_error->{$mz}->{'intensity'},"\t",$mass_error->{$mz}->{'error'},"\n";
	}
}

sub changeMH
{
	my $self=shift;
	my $H=$self->get_H_value();
	my $dta = $self->get_dta();
	my $dirname  = dirname($dta);
	my ($mz_hash,$charge,$mz)=@_;

	my $new_mz = ($mz-$H)*$charge+$H;

	if($new_mz>1500)
	{
		delete $mz_hash->{$mz};
		return 0;
	}

	if(!defined($mz_hash->{$mz}))
	{
		$mz_hash->{$mz}= 0;
	}
	$mz_hash->{$new_mz}=$mz_hash->{$mz};

	system(qq(cp $dirname/peaks/$mz $dirname/peaks/$new_mz));
	delete $mz_hash->{$mz};

}


sub deisotope_charge_1
{
	my ($self,$mz_hash,$selected_mz) = @_;
	my $parameter=$self->get_parameter();
	my $ppm = $parameter->{'charge12_ppm'};
	
	my ($peak_mz_low,$peak_mz_high) =  ($selected_mz - ($ppm/1000000)*$selected_mz, $selected_mz +($ppm/1000000)*$selected_mz);  

	foreach my $mz (sort {$a<=>$b} keys %$mz_hash)
	{
### check the next
		
		if($mz>$peak_mz_low && $mz<$peak_mz_high)
		{
			next if ($mz == $selected_mz);

			if(defined($mz_hash->{$mz}))
			{
				$mz_hash->{$selected_mz} += $mz_hash->{$mz};
### remove the un-deisotope peak				
				delete $mz_hash->{$mz};
			}

		}
	}

}	

1;