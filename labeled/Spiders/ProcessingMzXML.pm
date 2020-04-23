#!/usr/bin/perl

######### Simulation ##########################################
#                                                             #
#       **************************************************    #  
#       **** Coverter MzXML to dta files	          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::ProcessingMzXML;

use strict;
use warnings;
use File::Basename;
use Spiders::XMLParser;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT = qw(set_mzXML_file get_mzXML_file get_mzXMLfile_dirname set_dta_path get_dta_path set_converter get_converter set_log_file get_log_file mzXML2dta readmzXML get_last_scan_num get_total_scan_num get_index_array get_total_runtime set_parameter get_parameter get_ms1_rt generate_hash_dta generate_dta_file generate_MS2_file generate_MS1_file generate_MS1_file_unlabel generate_ms3_file strongest_prec strongest_prec_file strongest_prec_mz_array create_surveyhash);

sub new{
	my ($class,%arg)=@_;
    my $self = {
		_mzXML_file=>undef,
		_dta_path =>undef,
		_converter=>undef,
    };
    bless $self, $class;
	return $self;
}

sub set_mzXML_file
{
	my ($self,$mzXML_file)=@_;
	return $self->{_mzXML_file}=$mzXML_file;
}

sub get_mzXML_file
{
	my $self=shift;
	return $self->{_mzXML_file};
}

sub get_mzXMLfile_dirname
{
	my $self=shift;
	my $mzXMLfile = $self->get_mzXML_file();

	my @suffixlist=(".mzXML");
	my $mzXMLfile_dirname = dirname($mzXMLfile,@suffixlist);

	return $mzXMLfile_dirname;
}

sub set_dta_path
{
	my ($self,$dtapath) = @_;

	if(!(-e $dtapath))
	{
		system(qq(mkdir $dtapath >/dev/null 2>&1));
	}
	$self->{_dta_path}=$dtapath;
}

sub get_dta_path
{
	my ($self) = shift;
	return $self->{_dta_path};
}


sub set_converter
{
	my ($self,$converter)=@_;
	$self->{_converter}=$converter; 
}

sub get_converter
{
	my $self=shift;
	if(!defined($self->{_converter}))
	{
		$self->{_converter} = "/usr/local/bin/MzXML2Search ";	
	}
	return $self->{_converter};
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

sub mzXML2dta
{
	my $self=shift;
	my $mzXMLfile = $self->get_mzXML_file();
	my $dta_path = $self->get_dta_path();

	my $converter = $self->get_converter();

	print "  Converting the mzXML file to dta files\n";

	system(qq($converter $mzXMLfile >/dev/null 2>&1));

}

####### the following subroutine is used to process mzXML without relying on other program

sub readmzXML
{
	my ($self)=@_;
	my $mzXML = $self->get_mzXML_file();

	open (XML, "<$mzXML") || die "can not open the file";
	my $xml = new Spiders::XMLParser();
	my $indexOffset = $xml->get_IndexOffset(*XML); 
	my ($index_array, $last_scan) = $xml->get_IndexArray(*XML, $indexOffset);
	$self->{'_last_scan_num'} = $last_scan;
	$self->{'_total_scan_num'} = scalar (@$index_array);
	$self->{'_index_array'} = $index_array;
	
	my ($runtime) = $xml->get_RT(*XML, $$index_array[$last_scan-1]);
	$self->{'_total_runtime'} = $runtime;
	return *XML;
}

sub get_last_scan_num
{
	my $self=shift;
	return $self->{'_last_scan_num'};
}

sub get_total_scan_num
{
	my $self=shift;
	return $self->{'_total_scan_num'};
}

sub get_index_array
{
	my $self = shift;
	return $self->{'_index_array'};
}

sub get_total_runtime
{
	my $self=shift;
	return $self->{'_total_runtime'};
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

sub get_ms1_rt
{
	my ($self,$scan)=@_;
	my $mshash = $self->{'_mshash'};

	return $mshash->{'surveyhash'}->{$scan}->{'rt'};
}



sub generate_hash_dta
{
	my ($self, $ms_hash, $msms_hash, $mz_array, $param_hash) = @_;
	
	my $xml = new Spiders::XMLParser();
	my $parameter = $self->get_parameter();
	my $win = $parameter->{'isolation_window'};
	my $ms_hash4peakdection;
	my $LOG = $self->get_log_file();
	
	my %no_prec;

	my $XML = $self->readmzXML();
	my $index_array = $self->get_index_array();
	my $last_scan = $self->get_last_scan_num();

	my ($number, $scan_order, $prev_index, $prev_survey) = (0,0,0,0);
	my $found = 0;
	my $cycle = 0;

	my %MS1_peaks_for_search;
	print $LOG "  Generating MS1 and MS2 hashes          \n";
	
	foreach (@$index_array)
	{
		my $index = $_;
		my ($deleted);
		$number++;
		my $mslevel=0;


		if($parameter->{'first_scan_extraction'}>5)
		{
			next if ($number<($parameter->{'first_scan_extraction'}-5));
		}
		next if ($number>$parameter->{'last_scan_extraction'});			
		print "\r  Generating MS1 and MS2 hashes: $number of $last_scan scans          ";
		my ($rt) = $xml->get_RT(*XML, $index);

		($$msms_hash{$number}{'xml_mz'}, $$msms_hash{$number}{'xml_int'}, $$msms_hash{$number}{'xml_act'},$mslevel) = $xml->get_PrecursorMZINTACT(*XML, $index);
		$mslevel=$xml->get_MSLevel(*XML, $index);
        if($mslevel == 2)
        {
			next if ($number<$parameter->{'first_scan_extraction'});		

			next if ($prev_survey == 0);
			$$msms_hash{$number}{'survey'} = $prev_survey;

			my $mz = $$msms_hash{$number}{'xml_mz'};
			my $int = $$msms_hash{$number}{'xml_int'};			

			push(@{$MS1_peaks_for_search{$prev_index}},$mz);
			next if(!defined $mz);
			if ($mz == 0)
			{
				$$msms_hash{$number}{'unverified'} = 1;
				print "  The precursor ion for scan $number is not present in MS ...\n";
				print Dumper($$msms_hash{$number});
				print Dumper($$msms_hash{$number});
				exit;
			}
			$$msms_hash{$number}{'prec_int'} = $int;
			$$msms_hash{$number}{'prec_mz'} = sprintf("%.6f", $mz);
			$$msms_hash{$number}{'prec_MH'} = sprintf("%.6f", $mz);
			my @peaks_array;
			$xml->get_Peaks(*XML, \@peaks_array, $index);
			my $values = scalar(@peaks_array);
			my $peak_num = $values/2;
			my ($sum, $sum2) = (0,0);
			my @msms_mz;
			my @msms_int;	
			for (my $i=0;$i<$values;$i++)
			{
				
				$i++;
				push(@msms_mz,$peaks_array[$i-1]);
				push(@msms_int,$peaks_array[$i]);
				my $temp = $peaks_array[$i];
				$sum += $temp; 
				$sum2 += $temp**2;
			}
			$$msms_hash{$number}{'msms_mz'} = \@msms_mz;
			$$msms_hash{$number}{'msms_int'} = \@msms_int;
			
			my $total = $sum;
			$$msms_hash{$number}{'peak_num'} = $peak_num;
			$$msms_hash{$number}{'intensity'} = $total;
			$$msms_hash{$number}{'rt'} = $rt;
			
			my $scannum = $number;

		}
		elsif($mslevel == 1)
		{ # ms scan
		
		
		###### In the case of profile mode data, the MS1 is generated by 3D peak detection
		
			if($parameter->{'data_acquisition_mode'}==2)
			{
				$ms_hash->{'scan_order'} = $scan_order;
			}
			else
			{
				$self->create_surveyhash(*XML, $param_hash, \%{$$ms_hash{'surveyhash'}{$number}}, $index, $scan_order,\@{$$mz_array[$number]});
			}
			$$ms_hash{'orderhash'}{$scan_order} = $number;
			$$ms_hash{'surveyhash'}{$number}{'rt'} = $rt;
			$prev_survey = $number; $prev_index = $index;   $scan_order++;

#			$self->create_surveyhash(*XML, $param_hash, \%{$$ms_hash{$number}}, $index, $scan_order,\@{$$mz_array[$number]});
#			$$ms_hash{'orderhash'}{$scan_order} = $number;
#			$$ms_hash{$number}{'rt'} = $rt;
#			$prev_survey = $number; $prev_index = $index;   
			

			my @peakArray;	
			my $sizeArray = scalar(@peakArray);
			my @mz;
			my @int;
			for (my $i = 0; $i < $sizeArray; $i++) {
				$i++;
				push (@mz, $peakArray[$i - 1]);
				push (@int, $peakArray[$i]);
			}
			$ms_hash4peakdection->{$scan_order}{'RT'} = $rt;
			$ms_hash4peakdection->{$scan_order}{'scanNumber'} = $number;
			$ms_hash4peakdection->{$scan_order}{'mass'} = \@mz;
			$ms_hash4peakdection->{$scan_order}{'intensity'} = \@int;
			$scan_order++;			
		}
		elsif($mslevel == 3)
		{
			my @peaks_array;
			$xml->get_Peaks(*XML, \@peaks_array, $index);
			my $values = scalar(@peaks_array);

##################### generate the dta files ############################
			my $ms3_prec_mz = $$msms_hash{$number}{'xml_mz'};
			my $scannum = $number;
			$self->generate_ms3_file($scannum,$ms3_prec_mz,\@peaks_array);
		}
	}	

	$self->{'_mshash'} = $ms_hash;
	$self->{'_msmshash'} = $msms_hash;
	$self->{'_mshash4peakdetection'} = $ms_hash4peakdection;
}

sub generate_dta_file
{
	my ($self,$scannum,$mz,$peaks_array) = @_;
	my $H = 1.007276466812;
	my $dta_path = $self->get_dta_path();
	my $parameter = $self->get_parameter();	
	
	my $mzXML_filename = $self->get_mzXML_file();

	my @suffixlist=(".mzXML");
	my $Rawfile_basename = basename($mzXML_filename,@suffixlist);
	
	my $file = "$scannum" . ".MS2";

	my $values = scalar(@$peaks_array);

	return 0 if ($values<=10);
	open(FILE,">$dta_path/$file") || die "can not open the file\n";

	print FILE $mz," ",1,"\n";

	for (my $i=0;$i<$values;$i++)
	{
		$i++;			
		print FILE $$peaks_array[$i-1]," ",$$peaks_array[$i],"\n"; 
	}
	close(FILE);

}

sub generate_MS2_file
{
	my ($self,$msms_hash,$mz_array) = @_;
	my $H = 1.007276466812;
	my $dta_path = $self->get_dta_path();
	my $parameter = $self->get_parameter();	
	my $LOG = $self->get_log_file();
	
	my $MS1_peaks_sequenced;
	print "  Generating MS2 data files\n";
	print $LOG "  Generating MS2 data files\n";	
	foreach my $scannum (keys %$msms_hash)
	{
		next if ($scannum<$parameter->{'first_scan_extraction'});
		next if($$msms_hash{$scannum}{'xml_mz'}==0);
		
		my $charge = 1;
		my $file = "$scannum" . ".MS2";
		
		my $precursor_mz = $$msms_hash{$scannum}{'prec_mz'};
		my $MS1_scan = $$msms_hash{$scannum}{'survey'};

		
		next if(!defined($MS1_scan));
		next if(!defined($precursor_mz));
		
######		strongest_prec_file
		my $ms1_scan_file = "$dta_path/${MS1_scan}" . ".MS1";
		#my ($strongest_mz, $INT, $sumint) = $self->strongest_prec_file($ms1_scan_file,$precursor_mz);
		my ($strongest_mz, $INT) = $self->strongest_prec_mz_array($precursor_mz,\@{$mz_array->[$MS1_scan]});
		$$msms_hash{$scannum}{'prec_mz'} = $strongest_mz;
		
		my $msms_mz = $$msms_hash{$scannum}{'msms_mz'};
		my $msms_int = $$msms_hash{$scannum}{'msms_int'};		
		
		
		#my $mz=$$msms_hash{$scannum}{'prec_mz'};	
		#my $int=$$msms_hash{$scannum}{'prec_int'};	
		open(FILE,">$dta_path/$file") || die "can not open the file\n";

		#print FILE ($mz-$H)*1+$H," ",1,"\n";
		print FILE $strongest_mz," ",1,"\n";
		
		for (my $i=0;$i<=$#$msms_mz;$i++)
		{		
			print FILE $msms_mz->[$i]," ",$msms_int->[$i],"\n"; 
		}
		close(FILE);
#		push (@{$MS1_peaks_sequenced->{$MS1_scan}},$mz);
#		push (@{$MS1_peaks_sequenced->{$MS1_scan}},$int);		
		push (@{$MS1_peaks_sequenced->{$MS1_scan}},$strongest_mz);
		push (@{$MS1_peaks_sequenced->{$MS1_scan}},$INT);		

	}
	return $MS1_peaks_sequenced;
}

sub generate_MS1_file
{
	my ($self,$mz_array)=@_;
	my $dta_path = $self->get_dta_path();
	for(my $index=0;$index<=$#$mz_array;$index++) 
	{	
		next if($#{$mz_array->[$index]}<1);
		my $dtafile = "$dta_path/$index" . ".MS1";
		open (DTAFILE, ">", $dtafile);
		print DTAFILE  "1\t2\n";
		for (my $i = 0; $i <= $#{$mz_array->[$index]}; $i++) 
		{
			foreach my $mz (sort {$a<=>$b} keys %{$mz_array->[$index]->[$i]})
			{
				my $int = $mz_array->[$index]->[$i]->{$mz};
				print DTAFILE "$mz\t$int\n";
			}									
		}
		close(DTAFILE);

	}
}

sub generate_MS1_file_unlabel
{
	my ($self,$MS1_peaks_sequenced) = @_;
	my $dta_path = $self->get_dta_path();
	my @MS1_files = glob("$dta_path/*.MS1");

	my %sequenced_ms1;
	foreach my $scan (keys %$MS1_peaks_sequenced)
	{
		my $dtafile = "$dta_path/$scan" . ".MS1";

		for (my $i = 0; $i <= $#{$MS1_peaks_sequenced->{$scan}}; $i++) 
		{
			my ($mz,$int,$sum_int) = $self->strongest_prec_file($dtafile,$MS1_peaks_sequenced->{$scan}->[$i]);
			$sequenced_ms1{$dtafile}{$mz}=$int;	
			$i++;								
		}
	
	}
	
	foreach (@MS1_files)
	{
		if($sequenced_ms1{$_})
		{
			open (DTAFILE, ">", $_);
			print DTAFILE  "1\t2\n";

			foreach my $mz (sort {$a<=>$b} keys %{$sequenced_ms1{$_}})
			{
				next if($mz==0);
				print DTAFILE $mz," ",$sequenced_ms1{$_}{$mz},"\n";
			}
			close(DTAFILE);	
		}
		else
		{
			system(qq(rm -rf $_));
		}
	}	
	
}
sub generate_ms3_file
{
    my ($self, $scannum,$mz,$peaks_array) = @_;
	my $dta_path = $self->get_dta_path();
	my $Rawfile_basename = basename($dta_path);
	$Rawfile_basename =~ s/\..*//;

	my $file = "$Rawfile_basename\.$scannum\.$scannum" . ".ms3";
    my $values = scalar(@$peaks_array);


    open(FILE,">$dta_path/$file") || die "can not open the file\n";
    print FILE $mz," ",1,"\n";

    for (my $i=0;$i<$values;$i++)
    {
            $i++;
            print FILE $$peaks_array[$i-1]," ",$$peaks_array[$i],"\n";
    }
}

sub strongest_prec
{
	my ($self, $XML, $hash, $prev, $win) = @_;
	my $xml = new Spiders::XMLParser();
	my ($MZ, $INT) = (0, 0);

	# Find strongest precursor peak within +/- .003
	my $mz = $$hash{'xml_mz'};
	my ($low, $high) = ($mz-$win-.10, $mz+$win+.10); 

	my @peaks_array; 
	$xml->get_Peaks(*$XML, \@peaks_array, $prev);
	my $values = scalar(@peaks_array);
	my ($intensity, $sumint) = (0, 0);
	for (my $i=0;$i<$values;$i++)
	{
		my $mz = $peaks_array[$i];
		my $int = $peaks_array[$i+1];
		if ($mz < $low)
		{ 
			$i++;
			next;
		}
		last if ($mz > $high);
		$i++;
		$sumint+=$int;
		if ($int>$intensity)
		{
			$intensity = $int; 
			($MZ, $INT)  = ($peaks_array[$i-1], $intensity);	
		}
	}
	return ($MZ, $INT, $sumint); 
}

sub strongest_prec_file
{
	my ($self, $file, $mz) = @_;
	my $parameter = $self->get_parameter();	
	my $win = $parameter->{'isolation_window'};
	my ($MZ, $INT) = (0, 0);

	my ($low, $high) = ($mz-$win-.10, $mz+$win+.10); 
	my @peaks_array;
	return (0,0,0) if(!-e($file));
	open(INPUT,$file) || die "can not open the file: $!\n";
	while(<INPUT>)
	{
		chomp $_;
		my @data = split(/\t/,$_);
		push(@peaks_array,@data);
	}	
	
	my $values = scalar(@peaks_array);
	my ($intensity, $sumint) = (0, 0);
	for (my $i=0;$i<$values;$i++)
	{
		my $mz = $peaks_array[$i];
		my $int = $peaks_array[$i+1];
		if ($mz < $low)
		{ 
			$i++;
			next;
		}
		last if ($mz > $high);
		$i++;
		$sumint+=$int;
		if ($int>$intensity)
		{
			$intensity = $int; 
			($MZ, $INT)  = ($peaks_array[$i-1], $intensity);	
		}
	}
	return ($MZ, $INT, $sumint); 
}

sub strongest_prec_mz_array
{
	my ($self, $mz,$mz_array) = @_;
	my $parameter = $self->get_parameter();	
	my $win = $parameter->{'isolation_window'};
	my ($MZ, $INT) = (0, 0);

	my ($low, $high) = ($mz-$win-.10, $mz+$win+.10); 
	my @peaks_array;

	foreach my $mass (sort {$a<=>$b} keys %{$mz_array->[int($mz)]})
	{
		push (@peaks_array,$mass);
		push (@peaks_array,$mz_array->[int($mz)]->{$mass});		
	}
	
	my $values = scalar(@peaks_array);
	my ($intensity, $sumint) = (0, 0);
	for (my $i=0;$i<$values;$i++)
	{
		my $mz = $peaks_array[$i];
		my $int = $peaks_array[$i+1];
		if ($mz < $low)
		{ 
			$i++;
			next;
		}
		last if ($mz > $high);
		$i++;
		$sumint+=$int;
		if ($int>$intensity)
		{
			$intensity = $int; 
			($MZ, $INT)  = ($peaks_array[$i-1], $intensity);	
		}
	}
	return ($MZ,$INT); 
}

sub create_surveyhash{
	my ($self, $XML, $paramhash, $mshash, $index, $scan_order, $mz_array) = @_;
################### xusheng ###############################	

	my $parameter = $self->get_parameter();
	my $xml = new Spiders::XMLParser();
#	if ($parameter->{'data_acquisition_mode'} == 2){ return (0);}
	my $win = $parameter->{'isolation_window'};

#	my $win = 0.19; 
	if (!defined($win)){$win=0};
	# Get all peaks from survey scan and calculate background
	
	my ($sum, $sum2) = (0,0);	
	my @peaks_array;
	$xml->get_Peaks($XML, \@peaks_array, $index);
	my $values = scalar(@peaks_array);
	my $peak_num = $values/2;
	my %hash;	my @array;
	
	for (my $i=0;$i<$values;$i+=2)
	{
		my $origmz = sprintf("%.6f", $peaks_array[$i]);
		my $mz = $origmz; $mz =~ s/\..*//;
		my $int = $peaks_array[$i+1];
		$array[$mz]{$origmz} = $int;
		$hash{$origmz} = $int;
	}
	if ($win != 0)
	{
		my %temp = %hash;
		my $largest = 0;
		for my $origmz (sort {$hash{$b}<=>$hash{$a}} keys %hash)
		{
			next if (!defined($temp{$origmz}));
			my ($low, $high) = ($origmz-$win, $origmz+$win); 
			my ($lowint, $highint) = ($low, $high);
			for my $testmz (keys %{$array[$lowint]})
			{
				next if ($testmz > $high || $testmz < $low);
				next if ($testmz == $origmz);
				delete $temp{$testmz};
				#print "     removed $testmz\n";
			} 
			if ($lowint != $highint)
			{
				for my $testmz (keys %{$array[$highint]})
				{
					next if ($testmz > $high);
					next if ($testmz < $low || $testmz == $origmz);
					delete $temp{$testmz};
					#print "     removed $testmz\n";
				} 
			}
		}
		%hash = %temp;
	}
	while (my ($origmz, $int) = each %hash)
	{
		my $mz = $origmz;
		$mz =~ s/\..*//;
		$$mz_array[$mz]{$origmz} = $int;
	}
	$peak_num = scalar(keys %hash);

	# Delete scans without peaks
	if ($values <= 0){ return (0);}

	$$mshash{'scan_order'} = $scan_order; 

}

1;
