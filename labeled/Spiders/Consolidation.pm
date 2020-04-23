#!/usr/bin/perl

######### Simulation ##########################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################
# added the normalization function (v1.0.1) on 5/11/2012

package Spiders::Consolidation;
        
use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.01;
@ISA	 = qw(Exporter);
@EXPORT = qw(get_msms_mz get_msms_int get_keepnum set_parameter get_parameter set_H_value get_H_value get_strong_peak_window Consolidation Consolidation2 Normalize_Consolidation get_top_peak remove_neutral_losses read_dta_file write_dta_file check_pho_loss);

sub new{
	my ($class,%arg)=@_;
    my $self = {
        _msms_mz => $arg{'-msms_mz'},
        _msms_int => $arg{'-msms_int'},		
		_keepnum =>$arg{'-keepnum'},
    };
    bless $self, $class;
	return $self;
}

sub get_msms_mz
{
	my $self = shift;
	return $self->{_msms_mz};	
}

sub get_msms_int
{
	my $self = shift;
	return $self->{_msms_int};	
}

sub get_keepnum
{
	my $self = shift;
	return $self->{_keepnum};
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


sub get_strong_peak_window
{
	my $self=shift;
	return $self->{'_strong_peak_window'};
}

sub Consolidation
{
	my $self = shift;
	my @updated_msms_mz=();
	my $msms_mz = $self->get_msms_mz();

	my $msms_int = $self->get_msms_int();	
	
### generate a temp hash
	my %temp_hash=();
	for(my $i=0; $i<=$#$msms_mz; $i++)
	{
		$temp_hash{$msms_mz->[$i]}=$msms_int->[$i];
	}
	
	my $keep = $self->get_keepnum();

	my $startmz =  0;

	$startmz = (sort {$a <=> $b} @$msms_mz)[0];	

	$startmz = sprintf("%f", $startmz);
	my $nextmz = $startmz + 100; 
	my %temp;
	my $num = scalar(@$msms_mz);
	my $i = 0;
	my $j=0;
	my $win_num=0;
	my %keepers;
	foreach my $mz (sort {$a<=>$b} keys %temp_hash)
	{
###### $j is used to control the last data point	
		$j++;

		if($mz<=$nextmz)
		{
			$temp{$mz}=$temp_hash{$mz};
		}
	## goes to next window .. 	

		if($mz>$nextmz)
		{
			if((scalar keys %temp)>1)
			{
				my $k=0;
				%keepers=();
				for my $mz (reverse sort {$temp{$a}<=>$temp{$b}}  keys %temp)
				{
					$keepers{$mz} = $temp{$mz};
					$k++; 
					last if ($k==$keep);
				}
				for my $mz (sort {$a<=>$b} keys %keepers)
				{
					push(@updated_msms_mz,$mz);
				}
				$win_num++;
			}
			elsif((scalar keys %temp)==1)
			{
				for my $mz (keys %temp)
				{				
	
					push(@updated_msms_mz,$mz);
				}
				$win_num++;
			}			
			%temp=();
			$temp{$mz}=$temp_hash{$mz};			
			$nextmz += 100;				
		}
	}
	my $k=0;
	%keepers=();
	for my $mz (reverse sort {$temp{$a}<=>$temp{$b}}  keys %temp)
	{
		$keepers{$mz} = $temp{$mz};
		$k++; 
		last if ($k==$keep);
	}
	for my $mz (sort {$a<=>$b} keys %keepers)
	{
		push(@updated_msms_mz,$mz);
	}
	return \@updated_msms_mz;
}

sub Consolidation2
{
	my $self = shift;

	my $dta = $self->get_dta();
	my $keep = $self->get_keepnum();
	
}

sub Normalize_Consolidation
{
	my ($self)=@_;
	my $strong_peak_window = $self->get_strong_peak_window();
	my $dtafile = $self->get_dta();
	my $dtahash = $self->read_dta_file($dtafile);

	my $keep = $self->get_keepnum();
######## check if there are any peaks with pho neutral loss
	my $param = $self->get_parameter();
	my $pho_loss_num = 0;
	if($param->{'pho_neutral_loss'})
	{
		$pho_loss_num = $self->check_pho_loss($dtahash);
	}
##############################################################		
	my $top_peak = $self->get_top_peak($strong_peak_window);
	my $rank_strong_peak;
	my $j=0;
	foreach (reverse sort {$strong_peak_window->{$a}<=>$strong_peak_window->{$b}} keys %$strong_peak_window)
	{
		$rank_strong_peak->{$strong_peak_window->{$_}}=$j;
		$j++;
	}

	my $i=0;
	foreach my $mz (sort {$a<=>$b} keys %{$dtahash->{'ms2'}})
	{

		my $window_num = int($i/$keep);
		$i++;
		my $strong_peak_within_window = $strong_peak_window->{$window_num};
		next if (!defined($strong_peak_within_window));

		$dtahash->{'ms2'}->{$mz} = $dtahash->{'ms2'}->{$mz}*($top_peak/$strong_peak_within_window)*(1-0.01*$rank_strong_peak->{$strong_peak_within_window});

	}
	$self->write_dta_file($dtafile,$dtahash);
	return $pho_loss_num;
}

sub get_top_peak
{
	my ($self,$hash)=@_;
	my $max = 0;
	foreach my $value (keys %$hash ) {
		$max = $hash->{$value} if ($hash->{$value} > $max);
	}
	return $max;
}

sub remove_neutral_losses
{

}

sub read_dta_file
{
	my ($self,$dtafile) = @_;
	
	open (DTA, "$dtafile") || die "can not open the dta file: $dtafile\n";
	my %dtahash;
	my $line0 = <DTA>;
	my @data=split(/\s+/,$line0);

	$dtahash{'prec_mz'} = $data[0];
	$dtahash{'prec_charge'} = $data[1];

	while(<DTA>)
	{
		chomp $_;
		my ($mz,$int)=split(/\s+/,$_);
		$dtahash{'ms2'}{$mz}=$int;
	}
	close(DTA);
	return \%dtahash;
}

sub write_dta_file
{
	my ($self,$dtafile,$dtahash) = @_;

	if(-e $dtafile)
	{
		system(qq(rm $dtafile));
	}
	open (DTA, ">$dtafile") || die "can not open the dta file: $dtafile\n";
	print DTA $dtahash->{'prec_mz'}," ",$dtahash->{'prec_charge'},"\n";

	foreach my $mz (sort {$a<=>$b} keys %{$dtahash->{'ms2'}})
	{
		my $intensity = sprintf("%.6f",$dtahash->{'ms2'}->{$mz});
		$mz = sprintf("%.6f",$mz);
		print DTA $mz," ",$intensity,"\n";
	}
	close(DTA);
}

######### add on 6/10/2013 for checking the phosophorylation neutral loss 

sub check_pho_loss
{
	my ($self,$dtahash_orig) = @_;
	my $param = $self->get_parameter();
	my $pho_mass = $self->get_pho_neutral_loss();
	my $H = $self->get_H_value();
	
	my $prec_mz = $dtahash_orig->{'prec_mz'};
	my $charge = $dtahash_orig->{'prec_charge'};

	my $pho_peak = ($prec_mz - $H)	/ $charge + $H - $pho_mass / $charge;

	my $pho_loss_time = 0;

	foreach my $mz (keys %{$dtahash_orig->{'ms2'}})
	{
		next if (($mz-1)>$pho_peak);
		
		if(abs($pho_peak - $mz) < $param->{'frag_mass_tolerance'})
		{
			$pho_loss_time++;
		}
	}
	return $pho_loss_time;
}


 
1;