#!/usr/bin/perl

## Release date: 01/31/2015
## Release version: version 11.1.1
## Module name: Spiders::MassCorrection

##################################################################
##	Mass table calculated by Junmin Peng on 10/03/2014	##
##################################################################
## 	TMT reporters and y1-ions
##	 	126 = 126.1277259380
## 		127N = 127.1247608314
##		127C = 127.1310807758
##		128N = 128.1281156692
##		128C = 128.1344356136
##		129N = 129.1314705070
##		129C = 129.1377904514
##		130N = 130.1348253448
##		130C = 130.1411452892
##		131 = 131.1381801826
##		K = 376.2757362992
##		R = 175.1189521741
##	y1-ions of non-TMT
##		K = 147.1128041645
##		R = 175.1189521741
##	ESI MS1-ion
##		(Si(CH3)2O))6 + H+ = 445.1200245337

package Spiders::MassCorrection;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);


$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT = qw(set_log_file get_log_file massCorrection  partition  stdev  mean  getObservedMasses); 

sub new {
	my ($class,%arg)=@_;
    my $self = {};
    bless $self, $class;
	return $self;
}

sub set_log_file
{
	my ($self,$log)=@_;
	$self->{_log}=$log;	
}

sub get_log_file
{
	my ($self)=@_;
	return $self->{_log};	
}

sub massCorrection {
	##########################
	## Initialization	##
	##########################
	
	my ($self, $msHash, $ms2Hash, $mzArray, $correction_data, $params,$mass_shift_file) = @_;

	my $LOG = $self->get_log_file();
	my @options = split(/,/, $params->{'mass_correction'});
	my $option = shift (@options);
	my $manualMassShift = shift (@options);
	if (defined $manualMassShift) {
		$manualMassShift =~ s/\s//;
	}
	my $tolPpm = 50;				## Window-width for each reference ion
	my $thresholdPercentage = 0.05;			## Threshold ratio of the number of spectra used for the mass-shift calculation to the total number of spectra	 
							## e.g. $thresholdPercentage = 0.2 indicates that at least 20% of spectra (MS1 or MS2)
							## should be used to calculate mass shifts
	my $mass_shift_mean_all;									
	##################################
	## Mass-shift correction module	##
	##################################

	my %massShifts;
	my $nScans = scalar (@$mzArray);
	
	if ($option == 0) {
		##########################
		## No correction	##
		##########################
		print "  No mass-shift correction\n";
		print $LOG "  No mass-shift correction\n";		
		return ($ms2Hash, $mzArray);
	}
	elsif ($option == 1) 
	{
		##################################################
		## Mass-shift correction based on MS1 spectra	##
		##################################################
		foreach my $scanNumber (keys %$correction_data)	
		{
			if((scalar @{$correction_data->{$scanNumber}})>1)
			{
				my $min_diff = 1;
				my $select_mz=0;
				foreach (@{$correction_data->{$scanNumber}})
				{
					my $diff = $_ - 445.1200245337;
					if($min_diff > $diff)
					{
						$select_mz = $_;
					}
				}
				@{$correction_data->{$scanNumber}} = ($select_mz);				
			}
			$massShifts{$scanNumber} = ($correction_data->{$scanNumber}->[0]-445.1200245337)/445.1200245337 * 1e6;
		}
		
		##########################################################################################
		## Filter the calculated mass shifts							##
		## 1. More than 50% of spectra should be used to calculate mass shifts			##
		## 2. Remove 20% highest/lowest mass shifts to get a more reliable correction factor	## 
		##########################################################################################
		my $nMassShifts = scalar(keys %massShifts);
		my $defined_peaks = (scalar (keys %$msHash));
		my $measuredPercentage = $nMassShifts / (scalar (keys %$msHash));	
		## More than 5% of spectra should be used to calculate mass shifts
		## Otherwise, mass-shift correction will not be performed
		if ($measuredPercentage >= $thresholdPercentage) {
			printf("  %.2f%% of spectra (%d of %d) is used to calculate mass shifts\n", $measuredPercentage * 100, $nMassShifts, $defined_peaks);
			printf $LOG "  %.2f%% of spectra (%d of %d) is used to calculate mass shifts\n", $measuredPercentage * 100, $nMassShifts, $defined_peaks;
			## Filter the 20% of highest/lowest mass shifts
			my $numFiltered = int(0.1 * scalar(keys %massShifts));
			my $Up_Filtered = scalar(keys %massShifts) - $numFiltered;
			my $i=0;
			foreach my $key (sort {$massShifts{$a} <=> $massShifts{$b}} keys %massShifts)
			{
				if($i<$numFiltered)
				{
					delete $massShifts{$key};

				}
				if($i>$Up_Filtered)
				{
					delete $massShifts{$key};
				}
				$i++; 
			}

		} else {
			printf("  Only %.2f%% of spectra (%d of %d) is used to calculate mass shifts\n", 
					$measuredPercentage * 100, scalar(keys %massShifts), scalar (keys %$msHash));
			printf("  Mass correction will not be performed when less than %.2f%% of spectra is used\n", $thresholdPercentage * 100);
			print "  No mass-shift correction\n";
			print $LOG "  No mass-shift correction\n";
			my $total_scan_number = scalar keys %{$ms2Hash};
			for(my $i=0;$i<=$total_scan_number;$i++)
			{
				$mass_shift_mean_all->{$i}->{mass_shift}=0;
			}
			return ($msHash,$ms2Hash, $mzArray,$mass_shift_mean_all);
		
			print "  Please choose another option like,\n";
			print "    \"0\" (no correction) or\n";
			print "    \"3\" (manual correction with a specified mass-shift value)\n";
			print "  for the \"mass_correction\" parameter and run again\n";
			exit;
		}		
		
		##############################################################
		##  Use bins to correct the mass shift                   ##
		##############################################################
		my $bin_size = 100;

		if((scalar keys %massShifts)>1000)
		{
			$bin_size = int((scalar keys %massShifts)/10)+1;
		}

		my ($mass_shift_mean)=partition(\%massShifts,$bin_size);
		my $bin_num=scalar keys (%$mass_shift_mean);
		my $bin=0;
		print "  Use $bin_num bin to correct the mass shift\n";
		print $LOG "  Use $bin_num bin to correct the mass shift\n";
		
		foreach my $scan (sort {$a<=>$b} keys %$mass_shift_mean)
		{
			$bin++;
			printf ("  Mass-shift at bin $bin: mean = %.5f ppm\n", $mass_shift_mean->{$scan}->{mass_shift});
			printf  $LOG "  Mass-shift at bin $bin: mean = %.5f ppm\n", $mass_shift_mean->{$scan}->{mass_shift};			
			
		}	
=head		
		open(FILE,">$mass_shift_file");
		foreach my $scanNumber (sort {$a <=> $b} keys %massShifts)
		{
			print FILE $scanNumber,"\t",$massShifts{$scanNumber},"\n";
		}
		close(FILE);
=cut		
		my $last_scan = 0;
		foreach my $scan (sort {$a<=>$b} keys %$mass_shift_mean)
		{
			for(my $i=0;$i<=$scan;$i++)
			{
				next if(defined($mass_shift_mean_all->{$i}));
				$mass_shift_mean_all->{$i}->{mass_shift}=$mass_shift_mean->{$scan}->{mass_shift};

			}
			$last_scan = $scan;
		}
		#######
		# The last scan in mass_shift_mean hash is not the last scan mshash and mzArray, use the last_scan to fill the remain data
		####
		$nScans = scalar keys %{$ms2Hash} if($nScans<(scalar keys %{$ms2Hash}));
		if($last_scan<=$nScans)
		{
			for(my $i=$last_scan;$i<=$nScans;$i++)
			{
				$mass_shift_mean_all->{$i}->{mass_shift}=$mass_shift_mean->{$last_scan}->{mass_shift};
			}
		}
		##########################################################################
		## Mass-shift correction of all MS1 and MS2 spectra using the correction factor	##
		##########################################################################
		print  $LOG "  correcting MS1 scan\n";
		for (my $scanNumber = 0; $scanNumber < $nScans; $scanNumber++)
		{
			print "\r  correcting MS1 scan: $scanNumber";
			## Correction of MS1 masses
			if (defined $$mzArray[$scanNumber]) {				
				for (my $i = 0; $i < scalar(@{$$mzArray[$scanNumber]}); $i++) {
					next if (!defined $$mzArray[$scanNumber][$i]);

					## Keys of %{$$mzArray[$scanNumber][$i]} (correspond to mz values) need to be changed
					%{$$mzArray[$scanNumber][$i]} = map {
						my $newKey = $_ / (1 + $mass_shift_mean_all->{$scanNumber}->{mass_shift} / 1e6);
						$newKey => $$mzArray[$scanNumber][$i]{$_}
					} keys (%{$$mzArray[$scanNumber][$i]});
				}
			}
				
		}
		print "\n";
		print  $LOG "  correcting MS2 scan\n";
		foreach my $scanNumber (keys %{$ms2Hash})
		{
			print "\r  correcting MS2 scan: $scanNumber";		
			if (defined $$ms2Hash{$scanNumber}{'msms_mz'}) 
			{
				## Correction of MS2 masses
				if (defined $$ms2Hash{$scanNumber}{'prec_mz'}) {
					for (my $i = 0; $i < scalar(@{$$ms2Hash{$scanNumber}{'msms_mz'}}); $i++) {
						$$ms2Hash{$scanNumber}{'msms_mz'}[$i] = $$ms2Hash{$scanNumber}{'msms_mz'}[$i] / (1 + $mass_shift_mean_all->{$scanNumber}->{mass_shift} / 1e6);
					}
				}
			}
			if (defined $$msHash{$scanNumber}{'mass'}) {				
				for (my $i = 0; $i < scalar(@{$$msHash{$scanNumber}{'mass'}}); $i++)
				{
					$$msHash{$scanNumber}{'mass'}[$i] = $$msHash{$scanNumber}{'mass'}[$i] / (1 + $mass_shift_mean_all->{$scanNumber}->{mass_shift} / 1e6);

				}
			}			
        }
	
		print "\n  Mass-shift correction has been finished\n";
		print $LOG "\n  Mass-shift correction has been finished\n";		
		return ($msHash,$ms2Hash, $mzArray,$mass_shift_mean_all);
	}

}

##################
## Subroutines	##
##################

sub partition {
    my ($hash, $N) = @_; 
	my $i=0;
	my %hash_bin;
	my %hash_bin_mean;
	foreach my $scan (sort {$a<=>$b} keys %$hash)
	{
		push(@{$hash_bin{int($i/$N)}{'mass_shift'}},$hash->{$scan});
		$hash_bin{int($i/$N)}{'scan'}=$scan;
		$i++;
	}
	foreach my $key (keys %hash_bin)
	{
		 $hash_bin_mean{$hash_bin{$key}{'scan'}}{'mass_shift'}= mean(@{$hash_bin{$key}{'mass_shift'}});
	}
	return \%hash_bin_mean;
}

sub stdev {
	my (@x) = @_;
	my $n = scalar(@x);
	my $stdev;	
	if ($n > 1) {
		my $mean = mean(@x);
		my $sum = 0;
		foreach (@x) {
			$sum = $sum + ($_ - $mean) ** 2;
		}
		$stdev = sqrt($sum / ($n - 1));
	} else {
		$stdev = 0;
	}
	return ($stdev);
}

sub mean {
	my (@x) = @_;
	my $n = 0;
	my $sum = 0;
	foreach (@x) {
		next if(!defined($_));
		$sum = $sum + $_;
		$n++;
	}
	my $mean = $sum / $n;
	return ($mean);
}

sub getObservedMasses {
	my ($massRef, $intensityRef, $referenceMasses, $tolPpm) = @_;
	my @masses = @{$massRef};
	my @intensities = @{$intensityRef};
	my @obsMasses = ();
	my $nPeaks = scalar(@masses);
	my $nReferences = scalar(@{$referenceMasses});
	for (my $i = 0; $i < $nReferences; $i++) {
		my $obsMass = "";
		my $obsIntensity = 0;	## This parameter can control the selection of a low-intensity peak
		my $uL = $$referenceMasses[$i] + $$referenceMasses[$i] * $tolPpm / 1e6;
		my $lL = $$referenceMasses[$i] - $$referenceMasses[$i] * $tolPpm / 1e6;
		for (my $j = 0; $j < $nPeaks; $j++) {
			if ($masses[$j] >= $lL && $masses[$j] <= $uL) {
				if ($intensities[$j] > $obsIntensity) {
					$obsMass = $masses[$j];
					$obsIntensity = $intensities[$j];
				}
			}
			last if ($masses[$j] > $uL);
		}
		push (@obsMasses, $obsMass);
	}
	return (@obsMasses);
}

1;
