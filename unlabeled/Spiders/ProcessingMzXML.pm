#!/usr/bin/perl

package Spiders::ProcessingMzXML;
use strict;
use warnings;
use File::Basename;
use Spiders::XMLParser;
use vars qw($VERSION @ISA @EXPORT);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = qw(set_mzXML_file get_mzXML_file get_mzXMLfile_dirname set_dta_path get_dta_path set_converter get_converter set_log_file get_log_file mzXML2dta readmzXML get_last_scan_num get_total_scan_num get_index_array get_total_runtime set_parameter get_parameter get_ms1_rt generate_hash_dta generate_dta_file generate_MS2_file generate_MS1_file generate_MS1_file_unlabel generate_ms3_file strongest_prec strongest_prec_file strongest_prec_mz_array create_surveyhash);

sub new {
	my ($class, %arg)=@_;
	my $self = {
		_mzXMLFile => undef,
		_dtaPath => undef,
		_converter => undef,
    };
    bless $self, $class;
	return $self;
}

sub setDtaPath {
	my ($self, $dtaPath) = @_;
	if (!-e $dtaPath) {
		system (qq(mkdir $dtaPath >/dev/null 2>&1));
	}
	$self -> {_dtaPath} = $dtaPath;
}

sub getDtaPath {
	my ($self) = shift;
	return $self->{_dtaPath};
}

sub setLogFile {
	my ($self, $log) = @_;
	$self -> {_log} = $log;	
}

sub getLogFile {
	my ($self) = @_;
	return $self -> {_log};	
}

sub setMzXMLFile {
	my ($self, $mzXMLFile) = @_;
	return $self -> {_mzXMLFile} = $mzXMLFile;
}

sub getMzXMLFile {
	my $self = shift;
	return $self -> {_mzXMLFile};
}

sub setParameter {
	my ($self, $param) = @_;
	$self -> {'_parameter'} = $param;
}

sub getParameter {
	my $self = shift;
	return $self -> {'_parameter'};
}

sub generateDtaHash {
	my ($self, $ms1Hash, $ms2Hash, $mzArray, $paramHash) = @_;
	my $xml = new Spiders::XMLParser();
	my $parameter = $self -> getParameter();
	my $win = $parameter->{'isolation_window'};	
	my $LOG = $self -> getLogFile();
	my $XML = $self -> readMzXML();
	my $indexArray = $self -> getIndexArray();
	my $lastScan = $self -> getLastScanNum();
	
	my ($scanNum, $scanOrder, $prevIndex, $prevSurvey) = (0, 0, 0, 0);
	my $found = 0;
	my $cycle = 0;
	my $msLevel = 0;
	my $ms1ForPeakDetection;
	my %ms1PeaksForSearch;
	print $LOG "  Generating MS1 and MS2 hashes          \n";	
	foreach my $index (@$indexArray) {		
		$scanNum++;
		if ($parameter -> {'first_scan_extraction'} > 5) {
			next if ($scanNum < ($parameter -> {'first_scan_extraction'} - 5));
		}
		next if ($scanNum > $parameter -> {'last_scan_extraction'});			
		print "\r  Generating MS1 and MS2 hashes: $scanNum of $lastScan scans          ";
		my $rt = $xml -> getRT(*XML, $index);
		($$ms2Hash{$scanNum}{'xml_mz'}, $$ms2Hash{$scanNum}{'xml_int'}, $$ms2Hash{$scanNum}{'xml_act'}, $msLevel) = $xml -> getPrecursorMZINTACT(*XML, $index);
		$msLevel = $xml -> getMSLevel(*XML, $index);

		## Read individual spectra
        if ($msLevel == 2) {
			next if ($scanNum < $parameter -> {'first_scan_extraction'});
			next if ($prevSurvey == 0);
			$$ms2Hash{$scanNum}{'survey'} = $prevSurvey;
			my $mz = $$ms2Hash{$scanNum}{'xml_mz'};
			my $int = $$ms2Hash{$scanNum}{'xml_int'};			
			push (@{$ms1PeaksForSearch{$prevIndex}}, $mz);
			next if(!defined $mz);
			if ($mz == 0) {
				$$ms2Hash{$scanNum}{'unverified'} = 1;
				print "  The precursor ion for scan#$scanNum is not present in MS ...\n";
				print Dumper($$ms2Hash{$scanNum});
				print $LOG Dumper($$ms2Hash{$scanNum});
				exit;
			}
			$$ms2Hash{$scanNum}{'prec_int'} = $int;
			$$ms2Hash{$scanNum}{'prec_mz'} = sprintf("%.6f", $mz);
			$$ms2Hash{$scanNum}{'prec_MH'} = sprintf("%.6f", $mz);
			my @peaksArray;
			$xml -> getPeaks(*XML, \@peaksArray, $index);
			my $nPeaks = scalar(@peaksArray) / 2;
			my $sumIntensity = 0;
			my @ms2Mz;
			my @ms2Int;	
			for (my $i = 0; $i < scalar(@peaksArray); $i++) {
				$i++;
				push (@ms2Mz, $peaksArray[$i - 1]);
				push (@ms2Int, $peaksArray[$i]);
				$sumIntensity += $peaksArray[$i];
			}
			$$ms2Hash{$scanNum}{'msms_mz'} = \@ms2Mz;
			$$ms2Hash{$scanNum}{'msms_int'} = \@ms2Int;
			$$ms2Hash{$scanNum}{'peak_num'} = $nPeaks;
			$$ms2Hash{$scanNum}{'intensity'} = $sumIntensity;
			$$ms2Hash{$scanNum}{'rt'} = $rt;			
		} elsif ($msLevel == 1) {
			## In the case of profile mode data, the MS1 is generated by 3D peak detection
			if ($parameter -> {'data_acquisition_mode'} == 2) {
				$ms1Hash -> {'scan_order'} = $scanOrder;
			} else {
				$self -> createSurveyhash(*XML, $paramHash, \%{$$ms1Hash{'surveyhash'}{$scanNum}}, $index, $scanOrder,\@{$$mzArray[$scanNum]});
			}
			$$ms1Hash{'orderhash'}{$scanOrder} = $scanNum;
			$$ms1Hash{'surveyhash'}{$scanNum}{'rt'} = $rt;
			$prevSurvey = $scanNum; 
			$prevIndex = $index;
			$scanOrder++;





			## What does this step mean?
			## According to the following codes, @mz and @int will be empty vectors




			my @peakArray;	
			my $sizeArray = scalar(@peakArray);
			my @mz;
			my @int;
			for (my $i = 0; $i < $sizeArray; $i++) {
				$i++;
				push (@mz, $peakArray[$i - 1]);
				push (@int, $peakArray[$i]);
			}
			$ms1ForPeakDetection -> {$scanOrder}{'RT'} = $rt;
			$ms1ForPeakDetection -> {$scanOrder}{'scanNumber'} = $scanNum;
			$ms1ForPeakDetection -> {$scanOrder}{'mass'} = \@mz;
			$ms1ForPeakDetection -> {$scanOrder}{'intensity'} = \@int;
			$scanOrder++;			
		} elsif ($msLevel == 3) {
			my @peaksArray;
			$xml -> getPeaks(*XML, \@peaksArray, $index);
			my $values = scalar(@peaksArray);
			my $ms3PrecMz = $$ms2Hash{$scanNum}{'xml_mz'};			
			$self -> generateMs3File($scanNum,$ms3PrecMz,\@peaksArray);
		}
	}
	$self -> {'_msHash'} = $ms1Hash;
	$self -> {'_msmsHash'} = $ms2Hash;
	$self -> {'_msHash4peakdetection'} = $ms1ForPeakDetection;
}

sub readMzXML {
	my ($self) = @_;
	my $mzXML = $self -> getMzXMLFile();
	open (XML, "<", "$mzXML") || die "Cannot open the file";
	my $xml = new Spiders::XMLParser();
	my $indexOffset = $xml -> getIndexOffset(*XML); 
	my ($indexArray, $lastScan) = $xml -> getIndexArray(*XML, $indexOffset);
	$self -> {'_last_scan_num'} = $lastScan;
	$self -> {'_total_scan_num'} = scalar (@$indexArray);
	$self -> {'_index_array'} = $indexArray;	
	my ($rt) = $xml -> getRT(*XML, $$indexArray[$lastScan - 1]);
	$self -> {'_total_runtime'} = $rt;
	return *XML;
}

sub getIndexArray {
	my $self = shift;
	return $self -> {'_index_array'};
}

sub getLastScanNum {
	my $self = shift;
	return $self -> {'_last_scan_num'};
}

sub createSurveyhash {
	my ($self, $XML, $paramHash, $msHash, $index, $scanOrder, $mzArray) = @_;	
	my $parameter = $self -> getParameter();
	my $xml = new Spiders::XMLParser();
	my $win = $parameter -> {'isolation_window'};
	if (!defined($win)) {
		$win = 0
	};
	
	## Get all peaks from a survey scan and calculate the background
	my @peaksArray;
	$xml -> getPeaks($XML, \@peaksArray, $index);
	my %hash;
	my @array;	
	for (my $i = 0; $i < scalar(@peaksArray); $i += 2) {
		my $origMz = sprintf("%.6f", $peaksArray[$i]);
		my $mz = $origMz; 
		$mz =~ s/\..*//;
		my $int = $peaksArray[$i + 1];
		$array[$mz]{$origMz} = $int;
		$hash{$origMz} = $int;
	}
	if ($win != 0) {
		my %temp = %hash;
		my $largest = 0;
		for my $origMz (sort {$hash{$b} <=> $hash{$a}} keys %hash) {
			next if (!defined($temp{$origMz}));
			my ($low, $high) = ($origMz - $win, $origMz + $win); 
			my ($lowInt, $highInt) = ($low, $high);
			for my $testmz (keys %{$array[$lowInt]}) {
				next if ($testmz > $high || $testmz < $low);
				next if ($testmz == $origMz);
				delete $temp{$testmz};
			}
			if ($lowInt != $highInt) {
				for my $testmz (keys %{$array[$highInt]}) {
					next if ($testmz > $high);
					next if ($testmz < $low || $testmz == $origMz);
					delete $temp{$testmz};
				} 
			}
		}
		%hash = %temp;
	}	
	while (my ($origMz, $int) = each %hash) {
		my $mz = $origMz;
		$mz =~ s/\..*//;
		$$mzArray[$mz]{$origMz} = $int;
	}

	# Delete scans without peaks
	if (scalar(@peaksArray) <= 0) {
		return (0);
	}
	$$msHash{'scan_order'} = $scanOrder;
}

sub generateMS2File {
	my ($self, $ms2Hash, $mzArray) = @_;
	my $H = 1.007276466812;
	my $dtaPath = $self -> getDtaPath();
	my $parameter = $self -> getParameter();	
	my $LOG = $self -> getLogFile();
	
	my $ms1PeaksSequenced;
	print "  Generating MS2 data files\n";
	print $LOG "  Generating MS2 data files\n";	
	foreach my $scanNum (keys %$ms2Hash) {
		next if ($scanNum < $parameter -> {'first_scan_extraction'});
		next if ($$ms2Hash{$scanNum}{'xml_mz'} == 0);
		my $charge = 1;
		my $file = "$scanNum" . ".MS2";
		my $precursorMz = $$ms2Hash{$scanNum}{'prec_mz'};
		my $ms1Scan = $$ms2Hash{$scanNum}{'survey'};		
		next if (!defined($ms1Scan));
		next if (!defined($precursorMz));

		my ($strongestMz, $intensity) = $self -> strongestPrecMzArray($precursorMz, \@{$mzArray -> [$ms1Scan]});
		$$ms2Hash{$scanNum}{'prec_mz'} = $strongestMz;
		my $ms2MzArray = $$ms2Hash{$scanNum}{'msms_mz'};
		my $ms2IntArray = $$ms2Hash{$scanNum}{'msms_int'};
		open(FILE, ">", "$dtaPath/$file") || die "Cannot open the file\n";
		print FILE $strongestMz, " ", 1, "\n";
		for (my $i = 0; $i <= $#$ms2MzArray; $i++) {		
			print FILE $ms2MzArray -> [$i], " ", $ms2IntArray -> [$i], "\n"; 
		}
		close(FILE);
		push (@{$ms1PeaksSequenced->{$ms1Scan}}, $strongestMz);
		push (@{$ms1PeaksSequenced->{$ms1Scan}}, $intensity);
	}
	return $ms1PeaksSequenced;
}

sub generateMS1File {
	my ($self, $mzArray) = @_;
	my $dtaPath = $self -> getDtaPath();
	for (my $index = 0; $index < scalar(@$mzArray); $index++) {
		next if ($#{$mzArray -> [$index]} < 1);
		my $dtaFile = "$dtaPath/$index" . ".MS1";
		open (DTAFILE, ">", $dtaFile);
		print DTAFILE  "1\t2\n";
		for (my $i = 0; $i <= $#{$mzArray->[$index]}; $i++) {
			foreach my $mz (sort {$a <=> $b} keys %{$mzArray->[$index]->[$i]}) {
				my $int = $mzArray -> [$index] -> [$i] -> {$mz};
				print DTAFILE "$mz\t$int\n";
			}
		}
		close(DTAFILE);
	}
}


=head


sub get_mzXMLfile_dirname
{
	my $self=shift;
	my $mzXMLfile = $self->get_mzXML_file();

	my @suffixlist=(".mzXML");
	my $mzXMLfile_dirname = dirname($mzXMLfile,@suffixlist);

	return $mzXMLfile_dirname;
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


sub get_total_scan_num
{
	my $self=shift;
	return $self->{'_total_scan_num'};
}


sub get_total_runtime
{
	my $self=shift;
	return $self->{'_total_runtime'};
}





sub get_ms1_rt
{
	my ($self,$scan)=@_;
	my $msHash = $self->{'_msHash'};

	return $msHash->{'surveyhash'}->{$scan}->{'rt'};
}




sub generate_dta_file
{
	my ($self,$scanNum,$mz,$peaks_array) = @_;
	my $H = 1.007276466812;
	my $dta_path = $self->get_dta_path();
	my $parameter = $self->get_parameter();	
	
	my $mzXML_filename = $self->get_mzXML_file();

	my @suffixlist=(".mzXML");
	my $Rawfile_basename = basename($mzXML_filename,@suffixlist);
	
	my $file = "$scanNum" . ".MS2";

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
	my ($self,$ms2Hash,$mz_array) = @_;
	my $H = 1.007276466812;
	my $dta_path = $self->get_dta_path();
	my $parameter = $self->get_parameter();	
	my $LOG = $self->get_log_file();
	
	my $MS1_peaks_sequenced;
	print "  Generating MS2 data files\n";
	print $LOG "  Generating MS2 data files\n";	
	foreach my $scanNum (keys %$ms2Hash)
	{
		next if ($scanNum<$parameter->{'first_scan_extraction'});
		next if($$ms2Hash{$scanNum}{'xml_mz'}==0);
		
		my $charge = 1;
		my $file = "$scanNum" . ".MS2";
		
		my $precursorMz = $$ms2Hash{$scanNum}{'prec_mz'};
		my $MS1_scan = $$ms2Hash{$scanNum}{'survey'};

		
		next if(!defined($MS1_scan));
		next if(!defined($precursorMz));
		
######		strongest_prec_file
		my $ms1_scan_file = "$dta_path/${MS1_scan}" . ".MS1";
		#my ($strongestMz, $INT, $sumint) = $self->strongest_prec_file($ms1_scan_file,$precursorMz);
		my ($strongestMz, $INT) = $self->strongest_prec_mz_array($precursorMz,\@{$mz_array->[$MS1_scan]});
		$$ms2Hash{$scanNum}{'prec_mz'} = $strongestMz;
		
		my $msms_mz = $$ms2Hash{$scanNum}{'msms_mz'};
		my $msms_int = $$ms2Hash{$scanNum}{'msms_int'};		
		
		
		#my $mz=$$ms2Hash{$scanNum}{'prec_mz'};	
		#my $int=$$ms2Hash{$scanNum}{'prec_int'};	
		open(FILE,">$dta_path/$file") || die "can not open the file\n";

		#print FILE ($mz-$H)*1+$H," ",1,"\n";
		print FILE $strongestMz," ",1,"\n";
		
		for (my $i=0;$i<=$#$msms_mz;$i++)
		{		
			print FILE $msms_mz->[$i]," ",$msms_int->[$i],"\n"; 
		}
		close(FILE);
#		push (@{$MS1_peaks_sequenced->{$MS1_scan}},$mz);
#		push (@{$MS1_peaks_sequenced->{$MS1_scan}},$int);		
		push (@{$MS1_peaks_sequenced->{$MS1_scan}},$strongestMz);
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
    my ($self, $scanNum,$mz,$peaks_array) = @_;
	my $dta_path = $self->get_dta_path();
	my $Rawfile_basename = basename($dta_path);
	$Rawfile_basename =~ s/\..*//;

	my $file = "$Rawfile_basename\.$scanNum\.$scanNum" . ".ms3";
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



=cut

1;
