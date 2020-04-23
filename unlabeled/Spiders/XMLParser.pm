#!/usr/bin/perl

package Spiders::XMLParser;
use strict;
use warnings;
use MIME::Base64;
use Data::Dumper;

my %residMass =(A => 71.0788,
				B => 114.5962,
				C => 103.1388,
				D => 115.0886,
				E => 129.1155,
				F => 147.1176,
				G => 57.0519,
				H => 137.1411,
				I => 113.1594,
				K => 128.1741,
				L => 113.1594,
				M => 131.1926,
				N => 114.1038,
				O => 114.1472,
				P => 97.1167,
				Q => 128.1307,
				R => 156.1875,
				S => 87.0782,
				T => 101.1051,
				V => 99.1326,
				W => 186.2132,
				X => 113.1594,
				Y => 163.1760,
				Z => 128.6231);

sub new {
	my ($class) = @_;
	my $self = {};
	bless ($self, $class);
	return $self;
}

sub getRT {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $rttime;
	seek (XML, $scanIndex, 0);
	while (<XML>){
		next if (!/retentionTime=/);
		chomp;
		$rttime = $_;
		$rttime =~ s/.*retentionTime="PT([\d+\.]+)S".*/$1/o;
		last;
	}	
	return $rttime;
}

sub getIndexOffset {
	shift @_;
	(*XML) = @_;
	seek (XML, -120, 2);
	my ($indexOffset);
	while (<XML>){
		next if (!/<indexOffset>/);
		chomp;
		$indexOffset = $_;
		$indexOffset =~ s/\s+<indexOffset>(\d+)<\/indexOffset>.*/$1/o;
		last;
	}
	return $indexOffset;
}

sub getIndexArray {
	shift @_;
	(*XML, my $indexOffset) = @_;
	my @indexArray;
	my $lastscan = 0;
	seek (XML, $indexOffset, 0);
	while (<XML>){
		next if (/^\s+<index/);
		last if (/^\s+<\/index>.*/);
		chomp;
		next if (/scan/);
		$lastscan++;
		my $index = $_;
		$index =~ s/[\s\t]+\<offset id="\d+"[\s]*>(\d+)<\/offset>.*/$1/o;
		push (@indexArray, $index);
	}
	return (\@indexArray, $lastscan);
}

sub getEntireScanInfo {
	shift @_;
	(*XML, my $scanIndex) = @_;
	seek (XML, $scanIndex, 0);
	while(<XML>){
		last if (/\<\/scan\>/);
		print "$_\n";
	}
}

sub getPeaksinfo {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my ($nPeaks, $peakLine);
	seek (XML, $scanIndex, 0);
	while(<XML>) {
		next if (!/peaksCount/);
		chomp;
		$nPeaks = $_;
		$nPeaks =~ s/\s+peaksCount="(\d+)".*/$1/o;
		last;
	}
	while (<XML>) {
		if (/<peaks precision=\"32\">/ || /pairOrder=\"m\/z-int\">/) {
			chomp;
			$peakLine = $_;
			if (/<peaks precision="32">/) {
				$peakLine =~ s/\s+<peaks precision="32">([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			} else {
				$peakLine =~ s/\s+pairOrder="m\/z-int">([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			}
			last;
		} else {
			next;
		}
	}
	return ($nPeaks, $peakLine);
}

sub getMSLevel {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $msLevel;	
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		next if (!/msLevel=/);
		chomp;
		$msLevel = $_;
		$msLevel =~ s/\s+msLevel="(\d)".*/$1/o;
		last;
	}	
	return $msLevel;
}
sub getPrecursorMZINTACT {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $precMz; 
	my $precInt = 0; 
	my $precAct = "CID";
	my $msLevel = 2;                                                                                                                            
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		if ((/filterLine/)) {
			chomp;
			if ($_ =~ m/ms(.*)(\s\d+\.\d+)\@([a-z]+).*(\s\d+\.\d+)\@([a-z]+)/) {
				($msLevel,$precMz,$precAct) = ($1, $2, $3);
				$msLevel =~ s/\s+//g;
				$precMz =~ s/\s+//g;
				$precAct =~ tr/a-z/A-Z/;
				last;
			} elsif ($_ =~ m/ms(.*)(\s\d+\.\d+)\@([a-z]+)/) {
				($msLevel,$precMz,$precAct) = ($1, $2, $3);
				$msLevel =~ s/\s+//g;
				$precMz =~ s/\s+//g;
				$precAct =~ tr/a-z/A-Z/;
				last;
			} elsif($_ =~ / ms \[\d+\./) {
				$msLevel = 1;
				$precMz = 0;
				last;
			}
		}
	}
	return ($precMz, $precInt, $precAct, $msLevel);
}

sub getPrecursor {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $precMz;
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		next if (!/<precursorMz/);
		if (/<\/precursorMz>/) {
			chomp;
			$precMz = $_;
			if ($precMz =~ /activationMethod/) {
				$precMz =~ s/.*>(\d.+)<\/precursorMz>/$1/;			
			} elsif ($precMz =~ /precursorCharge/) {
				$precMz =~ s/\s+.*>([e\d\.\+]+)<\/precursorMz>.*/$1/;
			} else {
				$precMz =~ s/.*>(\d+)<\/precursorMz>.*/$1/;
			}
			last;
		} else {
			while (<XML>){
				chomp;
				exit if (!/collisionEnergy=/);
				$precMz = $_;
				$precMz =~ s/\s+collisionEnergy="[\d\.]+">([e\d\.\+]+)<\/precursorMz>.*/$1/;
				last;
			}
			last;
		}
	}                                                                                                                                                             
	return $precMz;
}

sub getPrecursorIntensity {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $precIntensity;
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		next if (!/<precursorMz/);
		chomp;
		$precIntensity = $_;
		$precIntensity =~ s/precursorIntensity="([e\d\.\-\+]+)".*/$1/;
		last;
	}
	return $precIntensity;
}

sub getPrecursorinfo {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $ce;
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		next if (!/collisionEnergy=/);
		chomp;
		$ce = $_;
		$ce =~ s/\s+collisionEnergy="(\d+)".*/$1/;
		last;
	}	
	my $precMz;
	while (<XML>) {
		next if (!/<precursorMz precursorIntensity=/);
		chomp;
		$precMz = $_;
		$precMz =~ s/\s+<precursorMz precursorIntensity="[e\d\.\-\+]+">([e\d\.\+]+)<\/precursorMz>.*/$1/;
		last;
	}	
	return ($ce, $precMz);
}

sub getBasePeakinfo {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my ($mz, $intensity, $rt);
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		if (/basePeakIntensity=/){
			chomp ($intensity = $_);
			$intensity =~ s/.*basePeakIntensity=\"([\d\.]+)\".*/$1/o;
			last;
		}
	    if (/retentionTime=/) {
			chomp($rt = $_);
			$rt =~ s/.*retentionTime="PT([\d\.]+)S".*/$1/o;
	    }
	    next if (!/basePeakMz=/);
	    chomp;
	    $mz = $_;
	    $mz =~ s/.*basePeakMz="([\d\.]+)".*/$1/o;
	}                                                                                                                                                             
	return ($mz, $intensity, $rt);
}

sub getBasePeakIntensity {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $intensity;
	seek (XML, $scanIndex, 0);
	while (<XML>){
		if (/basePeakIntensity=/) {
			chomp($intensity = $_);
			$intensity =~ s/.*basePeakIntensity=\"([e\d\.\+]+)\".*/$1/o;
			last;
		}
	}
	if ($intensity =~ /e/) {
		$intensity =~ s/([\d\.]+)e\+[0]+(\d+)/$1/;
		$intensity *= 10**$2;
	}
	return ($intensity);
}

sub getPeaksCount {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my ($peaksCount);
	seek (XML, $scanIndex, 0);
	while (<XML>){
		next if (!/peaksCount/);
		chomp;
		$peaksCount = $_;
		$peaksCount =~ s/\s+peaksCount="(\d+)".*/$1/o;
		last;
	}
	return $peaksCount;
}

sub getPeaksString {
	shift @_;
	(*XML, my $scanIndex) = @_;
	my $peakLine;                                                                                                                                          
	seek (XML, $scanIndex, 0);
	while (<XML>){
		next if (!/<peaks precision="32">/);
		chomp;
		$peakLine = $_;
		$peakLine =~ s/\s+<peaks precision="32">([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
		last;
	}
	return $peakLine;
}

sub getPeaks {
	shift @_;
	(*XML, my $peakArray, my $scanIndex) = @_;
	my ($peaks, $peakLine);
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		if (/<peaks\sprecision="32"/ || /pairOrder="m\/z-int"[\s]*>/) {
			chomp;
			next if ($_ =~ /<peaks precision="32"\Z/);
			$peakLine = $_;
			if (/<peaks precision="32">/) {
				$peakLine =~ s/\s+<peaks precision="32"[\s\w\W\d\=\"]+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			} else {
				$peakLine =~ s/\s+pairOrder="m\/z-int"[\s]*>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			}
			last;
		} elsif (/compressedLen.*>/) {
			chomp;
			$peakLine = $_;
			$peakLine =~ s/\s+compressedLen="[\d]"\s+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			last;
		} elsif (/contentType.*>/) {
			chomp;
			$peakLine = $_;
			$peakLine =~ s/\s+contentType=".*">([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			last;
		} else {
			next;
		}
	}
	## Base64 Decode
	$peaks = decode_base64($peakLine);
	my @hostOrder32 = unpack("N*", $peaks);
	for (@hostOrder32){
		my $float = unpack("f", pack("I", $_));
		push (@$peakArray, $float);
	}
}

sub getScanNum{
	shift @_;
	(*XML, my $scanIndex) = @_;
	seek (XML, $scanIndex, 0);
	my $scanNum;
	while (<XML>){
		next if (!/<scan\snum=\"\d+\"/);
		$_ =~ s/
//g;
		chomp;
		$scanNum = $_;
		$scanNum =~ s/<scan\snum=\"(\d+)\"//;
		$scanNum = $1;
		last;
	}
	return $scanNum;
}

=head
sub MS_surveytype{ #DMD Jan 28, 2008
	shift @_;
	(*XML, my $index) = @_;
	my ($smz, $emz) = getmzrange("", *XML, $index);
	my ($peaksum) = getPeaksCount("", *XML, $index);
	my $sourceval = $peaksum/($emz-$smz);
	return ($sourceval);
}

sub hybrid_check{
	shift @_;
	(*XML) = @_;
	my $hybrid = 0;
	my $indexOffset = getIndexOffset("",  *XML);
  my ($index_array, $last_scan) = getIndexArray("", *XML, $indexOffset);
	my ($smz, $emz, $smz2, $emz2);
	my ($peaksum, $peaksum2) = (0,0);
	#open (OUT, ">test_orbi.txt");
	#open (OUT2, ">test_ltq.txt");
	my $testnum = 0; my $round = 0;
	for (my $i=0; $i<$last_scan; $i++){
		my $index = $$index_array[$i];
		next if (getMSLevel("", *XML, $index) != 1);
		next if (!defined($$index_array[$i+1]));
		my $next = $$index_array[$i+1];	next if (getMSLevel("", *XML, $next) != 2);
		next if (!defined($$index_array[$i-1]));
		my $prev = $$index_array[$i-1]; return $hybrid if (getMSLevel("", *XML, $prev) == 2);
		$round++;
		my ($smz, $emz) = getmzrange("", *XML, $index); my ($peaksum) = getPeaksCount("", *XML, $index);
		my ($smz2, $emz2) = getmzrange("", *XML, $prev); my ($peaksum2) = getPeaksCount("", *XML, $prev);
		my $sourceval2 = ($peaksum)/($emz-$smz);
		my $sourceval1 = ($peaksum2)/($emz2-$smz2);
	#	print OUT2 "ltq: $index $smz $emz $peaksum $sourceval2\n";
#		print OUT "orbi: $prev $smz2 $emz2 $peaksum2 $sourceval1\n";
		$testnum += 1 if ($sourceval1 < .15 && $sourceval2 > .15);
		#last if $round == 5;
	}
	#print "$testnum $round\n";
	$testnum /= $round;
	$hybrid = sprintf("%.0f", $testnum);
	#print "$testnum\n";exit;
	return $hybrid;
	#exit;
}

sub MS_Source{ #DMD June 28, 2007
	shift @_;
	(*XML) = @_;
	my $indexOffset = getIndexOffset("",  *XML);
  my ($index_array, $last_scan) = getIndexArray("", *XML, $indexOffset);
	my ($smz, $emz);
	my ($num, $peaksum) = (0,0);
	my ($minpeak, $maxpeak) = (10000, 0);
	for (my $i=0; $i<$last_scan; $i++){
		my $index = $$index_array[$i];
		next if (getMSLevel("", *XML, $index) != 1);
		next if (!defined($$index_array[$i+1]));
		my $next = $$index_array[$i+1];
		next if (getMSLevel("", *XML, $next) != 2);
		my $scan = getscannum("", *XML, $index);
		if (!defined($smz)) {($smz, $emz) = getmzrange("", *XML, $index);}
		my ($peakcount) = getPeaksCount("", *XML, $index);
		$minpeak = $peakcount if ($peakcount<$minpeak);
		$maxpeak = $peakcount if ($peakcount>$maxpeak);
		$num++; $peaksum+=$peakcount;
	}
	my $sourceval = ($peaksum/$num)/($emz-$smz);
	#printf "Average # of peaks in Survey Scans = %d (scan range = %d)\n", $peaksum/$num, $emz-$smz;
	#printf "sum=$peaksum num=$num range=%d average=%d min=$minpeak max=$maxpeak final=$sourceval\n", $emz-$smz, $peaksum/$num;
	return ($sourceval);
}

sub getmzrange{ #DMD June 28, 2007
	shift @_;
	(*XML, my $scanIndex) = @_;
  my ($smz, $emz);

  seek (XML, $scanIndex, 0);
  while (<XML>){
    next if (!/startMz=/ && !/endMz=/);
    chomp;
		if (/startMz=/){
    	$smz = $_; $smz =~ s/.*startMz=\"([\d+\.]+)\".*/$1/o;
		} else {
    	$emz = $_; $emz =~ s/.*endMz=\"([\d+\.]+)\".*/$1/o;
    	last;
		}
  }

  return ($smz, $emz);
}

sub getruntime{
	shift @_;
	(*XML) = @_;
	seek (XML, 0, 0);
	
	my $runtime;
	while (<XML>){
		next if (!/endTime=/);
		$runtime = $_;
		$runtime =~ s/^\s+endTime=\"PT([\d\.]+)S\".*/$1/;
		last;
	}
	
	return $runtime;
}


sub getendpointrts{
	shift @_;
	(*XML) = @_;
	seek (XML, 0, 0);
	my ($start, $end);
	while (<XML>){
		if (/startTime/ || /endTime/){
			chomp;
			$_ =~ s/
//g;
			if (/startTime/){
				$start = $_;
				$start =~ s/^\s+startTime=\"PT([0-9\.]+)S\"/$1/;
				$start *= 60;
			} else {
				$end = $_;
				$end =~ s/^\s+endTime=\"PT([0-9\.]+)S\">/$1/;
				$end *= 60;
				return ($start, $end);
			}
		}
	}
}
=cut