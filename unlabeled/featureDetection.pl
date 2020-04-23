#!/usr/bin/perl
        
use strict;
use warnings;
use File::Basename;
use List::Util qw(min max);
use MIME::Base64;
#use Statistics::Lite qw(mean median mode stddev);
use Storable;

################
## Parameters ##
################
my ($paramFile, $inputFile, $subDir) = @ARGV;
my %params = getParams($paramFile);
my $isCentroided = 0;
if ($params{'data_acquisition_mode'} == 1) {
	$isCentroided = 1;
} elsif ($params{'data_acquisition_mode'} == 2) {
	$isCentroided = 0;
} else {
	die ("Please set the proper 'data_acquisition_mode' parameter");
}
my $gap = $params{'skipping_scans'};	## Number of skipped scans (>= 0)
my $matchPpm = $params{'mass_tolerance_peak_matching'};	## The tolerance of mass difference (ppm) when growing 3D peaks
my $firstScanExtraction = $params{'first_scan_extraction'};
my $lastScanExtraction = $params{'last_scan_extraction'};
my $snRatio = $params{'signal_noise_ratio'};
my $minPeakIntensity = $params{'min_peak_intensity'};

my $scanWindow = $gap + 1;	## Number of neighboring MS1 scans for 3D peaks (= gap + 1)
my $nCentroidPoints = 3;	## The number of raw data points for calculating centroid m/z (in case of profile mode data)
my $maxPctRtRange = 100;
	
####################
## Get MS spectra ##
####################
my %ms = getMsHash($inputFile, $firstScanExtraction, $lastScanExtraction);
print "  Generation of 3D-peaks (Be patient, it takes time)\n";

###############################
## Mass ranges for each scan ##
###############################
my %massRanges;
my @minMass;
my @maxMass;
my $msCount = scalar(keys %ms);
for (my $i = 0; $i < $msCount; $i++) {
	push (@minMass, min(@{$ms{$i}{'mass'}}));
	push (@maxMass, max(@{$ms{$i}{'mass'}}));
}
%massRanges = ("minMass" => \@minMass, "maxMass" => \@maxMass);

###########################
## MaxQuant main routine ##
###########################
my @features;
my $nFeatures = -1;
my %noisyScan;
my @cache;
my $oldMinInd = -1;
my $oldMaxInd = -1;
if ($gap < 0) {
	die "\"scanWindow\" parameter should be greater than 0\n";
}
for (my $i = 0; $i < $msCount; $i++) {
	## Put former and latter scans to cache (temporary) array 
	my $minInd = max(0, $i - $gap - 1);
	my $maxInd = min($msCount - 1, $i + $gap + 1);
	if ($i == 0) {
		for (my $j = 0; $j <= $maxInd; $j++) {
			my ($msCentroidList,$intensityThreshold) = detectPeaks($nCentroidPoints, $ms{$j}, $isCentroided);
			$$msCentroidList{'scanNumber'} = $ms{$j}{'scanNumber'};
			$$msCentroidList{'scanIndex'} = $j;
			$$msCentroidList{'RT'} = $ms{$j}{'RT'};
			push @cache, $msCentroidList;
		}
	} else {
		for (my $j = $oldMinInd; $j < $minInd; $j++) {
			shift @cache;
		}
		for (my $j = $oldMaxInd + 1; $j <= $maxInd; $j++) {
			my ($msCentroidList,$intensityThreshold) = detectPeaks($nCentroidPoints, $ms{$j}, $isCentroided);
			$$msCentroidList{'scanNumber'} = $ms{$j}{'scanNumber'};
			$$msCentroidList{'scanIndex'} = $j;
			$$msCentroidList{'RT'} = $ms{$j}{'RT'};
			push @cache, $msCentroidList;
		}
	}		

	#################################################################################
	## Reduction step                                                              ##
	## For each scan, take the MS centroid list and select (i.e. reduce) the peaks ##
	## which can form a 3D-peak with other peaks in neighboring scans              ##
	## Reduction is done in forward and backward direction                         ## 
	#################################################################################
	my %p = %{$cache[$i - $minInd]};	## A MS centroid list under consideration
	my $pCount = @{$p{'peakCenterMass'}};
	my @valids;
	my $count = 0;
	for (my $j = 0; $j < $pCount; $j++) {
		my $cm = $p{'peakCenterMass'}[$j];
		my $match = 0;
		my $ntries = 0;
		for (my $k = $i - 1; $k >= $minInd; $k--) {
			if ($cm >= ($massRanges{'minMass'}[$k]) && $cm <= ($massRanges{'maxMass'}[$k])) {
				my %q = %{$cache[$k - $minInd]};
				my $ind;
				($match, $ind) = getClosestValue(\%q, $cm, $matchPpm);
				if ($match) {
					last;
				}
				$ntries++;
				if ($ntries > $scanWindow) {
					last;
				}
			}
		}
		if (!$match) {
			$ntries = 0;
			for (my $k = $i + 1; $k <= $maxInd; $k++) {
				next if (!defined($massRanges{'minMass'}[$k]));	
				next if (!defined($massRanges{'maxMass'}[$k]));	
				if ($cm >= ($massRanges{'minMass'}[$k]) && $cm <= ($massRanges{'maxMass'}[$k])) {
					my %q = %{$cache[$k - $minInd]};
					my $ind;
					($match, $ind) = getClosestValue(\%q, $cm, $matchPpm);
					if ($match) {
						last;
					}
					$ntries++;
					if ($ntries > $scanWindow) {
						last;
					}
				}
			}
		}
		if ($match) {
			push (@valids, $j);
		}
	}
	my ($reducedRef, $noiseRef) = extractMsCentroidList(\%p, \@valids);
	my %reduced = %$reducedRef;
	my %noise = %$noiseRef;
	my @noiseArray = @{$noise{peakIntensity}};
	my @noiseUpdated;
	foreach my $n_value (@noiseArray) {
		if (defined($n_value) and $n_value ne "") {
			push (@noiseUpdated, $n_value);
		}
	}
	my @noiseSort = sort {$b <=> $a} @noiseUpdated;
	splice(@noiseSort, 0, int(0.5 * $#noiseSort));
	$noisyScan{$noise{scanNumber}} = median(@noiseSort);
	
	#########################################################################################
	## Generation of 3D peaks by combining peaks within the specified mass tolerance range ##
	## in the former scans (+/- scanWindow)                                                ##
	## 3D-peak generation is performed only in backward direction                          ##
	#########################################################################################	
	$cache[$i - $minInd] = \%reduced;
	my $reducedCount = @{$reduced{'peakCenterMass'}};
	for (my $j = 0; $j < $reducedCount; $j++) {
		my $cm = $reduced{'peakCenterMass'}[$j];
		my $scan = $reduced{'scanNumber'};
		my $Intensity = $reduced{'peakIntensity'}[$j];
	}
		
	for (my $j = 0; $j < $reducedCount; $j++) {
		my $cm = $reduced{'peakCenterMass'}[$j];
		my $match = 0;
		my $ntries = 0;
		my @matchedPeaksInds;
		for (my $k = $i - 1; $k >= $minInd; $k--) {	## Search for previous scans stored in a cache array
			if ($cm >= ($massRanges{'minMass'}[$k]) && $cm <= ($massRanges{'maxMass'}[$k])) {
				my %q = %{$cache[$k - $minInd]};
				my ($matchIndicator, $ind) = getClosestValue(\%q, $cm, $matchPpm);
				## $matchIndicator = 1 means that the j-th (reduced) peak in the i-th scan
				## can form a 3D-peak with $ind-th (reduced) peak in the previous scan (%q)
				if ($matchIndicator) {
					push (@matchedPeaksInds, @{$q{'featureIndex'}}[@{$ind}]);
					$match = 1;
				}
			}
		}
		if ($match) {
			@matchedPeaksInds = unique(@matchedPeaksInds);
			my $featureInd;
			if ($#matchedPeaksInds > 0) {	## Match to the peaks in mutliple previous scans
				$featureInd = min(@matchedPeaksInds);
				foreach my $matchedPeaksInd (@matchedPeaksInds) {
					# Processing of singleton features
					if ($matchedPeaksInd ne $featureInd) {
						my $peakTobeRemoved = $features[$matchedPeaksInd];
						push (@{$features[$featureInd]{'centerMz'}}, @{$$peakTobeRemoved{'centerMz'}});
						push (@{$features[$featureInd]{'intensity'}}, @{$$peakTobeRemoved{'intensity'}});
						push (@{$features[$featureInd]{'scanIndex'}}, @{$$peakTobeRemoved{'scanIndex'}});
						push (@{$features[$featureInd]{'scanNumber'}}, @{$$peakTobeRemoved{'scanNumber'}});
						push (@{$features[$featureInd]{'scanRT'}}, @{$$peakTobeRemoved{'scanRT'}});
					
						# Revise the information in cache
						foreach my $scanInd (@{$features[$matchedPeaksInd]{'scanIndex'}}) {
							for (my $l = 0; $l <=$#cache; $l++) {
								if ($cache[$l]{'scanIndex'} == $scanInd) {
									for (my $m = 0; $m <= $#{$cache[$l]{'featureIndex'}}; $m++) {
										if ($cache[$l]{'featureIndex'}[$m] == $matchedPeaksInd) {
											$cache[$l]{'featureIndex'}[$m] = $featureInd;
										}
									}
								}
							}
						}
						delete $features[$matchedPeaksInd];
					}
				}
			} else {
				$featureInd = shift @matchedPeaksInds;
			}
			push (@{$cache[$i - $minInd]{'featureIndex'}}, $featureInd);	# dummy value for keeping array size
			push (@{$features[$featureInd]{'centerMz'}}, $reduced{'peakCenterMass'}[$j]);
			push (@{$features[$featureInd]{'minMz'}}, $reduced{'peakMinMass'}[$j]);
			push (@{$features[$featureInd]{'maxMz'}}, $reduced{'peakMaxMass'}[$j]);
			push (@{$features[$featureInd]{'intensity'}}, $reduced{'peakIntensity'}[$j]);
			push (@{$features[$featureInd]{'scanIndex'}}, $reduced{'scanIndex'});
			push (@{$features[$featureInd]{'scanNumber'}}, $reduced{'scanNumber'});
			push (@{$features[$featureInd]{'scanRT'}}, $reduced{'RT'});
		}
		if (!$match) {
			if ($i < $msCount) {
				$nFeatures++;
				push (@{$cache[$i - $minInd]{'featureIndex'}}, $nFeatures);
				push (@{$features[$nFeatures]{'centerMz'}}, $reduced{'peakCenterMass'}[$j]);
				push (@{$features[$nFeatures]{'minMz'}}, $reduced{'peakMinMass'}[$j]);
				push (@{$features[$nFeatures]{'maxMz'}}, $reduced{'peakMaxMass'}[$j]);
				push (@{$features[$nFeatures]{'intensity'}}, $reduced{'peakIntensity'}[$j]);
				push (@{$features[$nFeatures]{'scanIndex'}}, $reduced{'scanIndex'});
				push (@{$features[$nFeatures]{'scanNumber'}}, $reduced{'scanNumber'});
				push (@{$features[$nFeatures]{'scanRT'}}, $reduced{'RT'});
			}
		}
	}
	$oldMinInd = $minInd;
	$oldMaxInd = $maxInd;
	my $currentCount = $i + 1;
	print "  Generation of 3D-peaks: traversing scan#$currentCount out of $msCount\r";	
}

## Peak re-numbering (because some peaks (i.e. singleton 3D peaks) are 'undef')
my @tempFeatures;
for (my $i = 0; $i < scalar(@features); $i++) {
	if (defined $features[$i]) {
		push (@tempFeatures, $features[$i]);
	}		
}
@features = @tempFeatures;
undef @tempFeatures;

##############################################################
## Filtering step                                           ##
## A 3D-peaks may contain multiple peaks from one scan      ##
## In this case, choose only one with the largest intensity ##
##############################################################	
for (my $i = 0; $i < scalar(@features); $i++) {
	if (scalar(@{$features[$i]{'scanNumber'}}) != scalar(unique(@{$features[$i]{'scanNumber'}}))) {		
		my %tempHash;
		for (my $j = 0; $j < scalar(@{$features[$i]{'scanNumber'}}); $j++) {
			if (defined $tempHash{$features[$i]{'scanNumber'}[$j]}) {
				my $currentIntensity = $features[$i]{'intensity'}[$j];
				if ($currentIntensity > $tempHash{$features[$i]{'scanNumber'}[$j]}{'intensity'}) {
					$tempHash{$features[$i]{'scanNumber'}[$j]}{'intensity'} = $currentIntensity;
					$tempHash{$features[$i]{'scanNumber'}[$j]}{'index'} = $j;
				}
			} else {
				$tempHash{$features[$i]{'scanNumber'}[$j]}{'intensity'} = $features[$i]{'intensity'}[$j];
				$tempHash{$features[$i]{'scanNumber'}[$j]}{'index'} = $j;
			}
		}
		my @uniqueIndex;
		foreach my $key (sort {$a <=> $b} keys %tempHash) {
			push (@uniqueIndex, $tempHash{$key}{'index'});
		}			
		@{$features[$i]{'centerMz'}} = @{$features[$i]{'centerMz'}}[@uniqueIndex];			
		@{$features[$i]{'intensity'}} = @{$features[$i]{'intensity'}}[@uniqueIndex];
		@{$features[$i]{'maxMz'}} = @{$features[$i]{'maxMz'}}[@uniqueIndex];
		@{$features[$i]{'minMz'}} = @{$features[$i]{'minMz'}}[@uniqueIndex];
		@{$features[$i]{'scanIndex'}} = @{$features[$i]{'scanIndex'}}[@uniqueIndex];
		@{$features[$i]{'scanNumber'}} = @{$features[$i]{'scanNumber'}}[@uniqueIndex];
		@{$features[$i]{'scanRT'}} = @{$features[$i]{'scanRT'}}[@uniqueIndex];
	}
}

############################################################
## Generation of feature files (temporary and final ones) ##
############################################################
print "  Generating feature tables (it may take a while)\n";
my $tmpFeatureFile = makeTemporaryFeatureTable($inputFile, \@features, \%noisyScan, $subDir, \%params);
my $featureFile = makeFeatureTable($inputFile, $tmpFeatureFile);
print "  Feature file has been generated for $inputFile\n";

#################
## Subroutines ##
#################
sub getParams {
	my ($file) = @_;
	if (!-e $file) {
		die "Cannot open a parameter file. Please check the path of the parameter file\n";
	}
	my %params;
	my $nInputFiles = 0;
	open (PARAMS, "<", $file);
	while (<PARAMS>) {
		next if (!/\=/ || /^\#/ || /^\s+\#/);
		$_ =~ s/\s+//g;
		$_ =~ s/\#.*//;
		my ($key, $value)  = split("\=", $_);
		if ($key eq "input_file") {
			$params{$key}[$nInputFiles] = $value;
			$nInputFiles++;
		} else {
			$params{$key} = $value;
		}
				
	}
	close (PARAMS);
	return (%params);
}

sub unique {
	my %seen;
	my @uniq =  grep { !$seen{$_}++ } @_;
	return (@uniq);
}

sub calcLocalMinPositions {
	my @y = @_; ## smoothed intensity
	my @minPos;
	my $length = @y;
	for (my $i = 2; $i < $length - 2; $i++) {
		my $b2 = $y[$i - 2];
		my $b1 = $y[$i - 1];
		my $x = $y[$i];
		my $a1 = $y[$i + 1];
		my $a2 = $y[$i + 2];
		if (isMin($b2, $b1, $x, $a1, $a2)) {
			push (@minPos, $i);
		}
	}
	my @minY = @y[@minPos];
	my @yOrderedIndex = sort {$minY[$a] <=> $minY[$b]} 0..$#minY;
	@minPos = @minPos[@yOrderedIndex];
	return (@minPos);
}

sub isMin {
	my ($b2, $b1, $x, $a1, $a2) = @_;
	my $isMin = 0;
	if ($x < $b1 && $x < $a1) {
		$isMin = 1;
	} elsif ($x == $b1 && $x < $b2 && $x < $a1) {
		$isMin = 1;
	} elsif ($x == $a1 && $x < $a2 && $x < $b1) {
		$isMin = 1;
	} elsif ($x < $b2 && $x == $b1 && $x == $a1 && $x < $a2) {
		$isMin = 1;
	} else {
		$isMin = 0;
	}
	return $isMin;	
}

sub extractMsCentroidList {
	my ($hash, $index) = @_;
	my %reduced;
	my $subInd = 0;
	my %noise;
	foreach my $key (keys %{$hash}) {
		if ($key eq "scanIndex" || $key eq "scanNumber" || $key eq "RT") {
			$reduced{$key} = $$hash{$key};
			$noise{$key} = $$hash{$key};
		} elsif ($key eq "peakCenterMassHash") {
			next;
		} elsif ($key eq "peakCenterMassIndexHash") {
			next;
		} elsif ($key eq "peakCenterMass") {
			my @sub = @{$$hash{$key}}[@$index];
			my @noise_array = @{$$hash{$key}};
			delete @noise_array[@$index];
			$reduced{$key} = \@sub;
			$noise{$key} = \@noise_array;
			$reduced{'peakCenterMassHash'} = {};
			$reduced{'peakCenterMassIndexHash'} = {};
			foreach my $subMass (@sub) {
				push (@{$reduced{'peakCenterMassHash'}{int($subMass)}}, $subMass);
				$reduced{'peakCenterMassIndexHash'}{$subMass} = $subInd;
				$subInd++;
			}
		} else {
			my @sub = @{$$hash{$key}}[@$index];
			$reduced{$key} = \@sub;
			my @noiseArray = @{$$hash{$key}};
			delete @noiseArray[@$index];
			$noise{$key} = \@noiseArray;			
		}
	}
	return (\%reduced, \%noise);
}

sub getClosestValue {
	my ($hash, $mass, $tol) = @_;
	my $lL = $mass - ($mass * $tol / 1e6);
	my $uL = $mass + ($mass * $tol / 1e6);
	my @searchKeys = (int($lL)..int($uL));
	
	my @returnIndex;
	my $intensity = 0;
	my $match = 0;
	
	foreach my $key (@searchKeys) {
		next if (!defined $$hash{'peakCenterMassHash'}{$key});
		for (my $i = 0; $i < scalar(@{$$hash{'peakCenterMassHash'}{$key}}); $i++) {
			my $value = $$hash{'peakCenterMassHash'}{$key}[$i];			
			if ($value >= $lL && $value <= $uL) {
				## returnIndex can be an array (there might be multiple indices in one scan)
				push (@returnIndex, $$hash{'peakCenterMassIndexHash'}{$value});
				$match = 1;
			}
		}
	}
	return ($match, \@returnIndex);
}

sub detectPeaks {
	## Input is a single MS spectrum (mz values and corresponding intensities)
	my ($numCentroidPoints, $msHash, $isCentroided) = @_;
	my @peakIntensity;
	my @peakCenterMass;
	my @peakMinMass;
	my @peakMaxMass;
	my @peakIndex;
	my %peakCenterMassHash;
	my %peakCenterMassIndexHash;
	my $peakCount = 0;
	my $intensityThreshold = 0;
	my $counts = @{$$msHash{'intensity'}};
	if ($isCentroided == 0) {	## Profile mode
		for (my $i = 2; $i < $counts - 2; $i++) {
			if ($$msHash{'intensity'}[$i] > 0) {
				# Consider 2 points before and after the point of interest (x), i.e. 5 point-window
				my $b2 = $$msHash{'intensity'}[$i - 2];
				my $b1 = $$msHash{'intensity'}[$i - 1];
				my $x = $$msHash{'intensity'}[$i];
				my $a1 = $$msHash{'intensity'}[$i + 1];
				my $a2 = $$msHash{'intensity'}[$i + 2];
				if ($x >= $intensityThreshold) {
					if (isMax($b2, $b1, $x, $a1, $a2)) {
						# If x is the local maximum within a 5-point window, lower and upper bounds for a peak will be explored
						# Refer Figure 1a and b in the paper (Cox and Mann, Nature Biotechnology, 2008, 26: 1367-72)
						my $minInd = calcMinPeakIndex($i, $$msHash{'intensity'});	# Search for a lower bound of the peak
						my $maxInd = calcMaxPeakIndex($i, $$msHash{'intensity'});	# Search for a upper bound of the peak
						if (($maxInd - $minInd) > 2) {
							my ($calculatedIntensity, $calculatedMass) = calcCenterMass($numCentroidPoints, $minInd, $i, $maxInd, $$msHash{'mass'}, $$msHash{'intensity'});
							if ($calculatedMass =~ /^?\d+\.?\d*$/) {
								$peakCenterMass[$peakCount] = $calculatedMass;
								$peakIntensity[$peakCount] = $calculatedIntensity;
								push (@{$peakCenterMassHash{int($calculatedMass)}}, $calculatedMass);
								$peakCenterMassIndexHash{$calculatedMass} = $peakCount;
								$peakIndex[$peakCount] = [$minInd..$maxInd];
								$peakCount++;
							}
						}
					}
				}
			}
		}		
	} else {	## Centroid mode
		@peakCenterMass = @{$$msHash{'mass'}};
		@peakIntensity = @{$$msHash{'intensity'}};
		for (my $i = 0; $i < scalar(@peakCenterMass); $i++) {
			$peakIndex[$i] = [$i + 1];
			$peakCenterMassIndexHash{$peakCenterMass[$i]} = ($i + 1);
			push (@{$peakCenterMassHash{int($peakCenterMass[$i])}}, $peakCenterMass[$i]);
		}
	}
	my %msCentroidList = ('peakCenterMass' => \@peakCenterMass, 
						'peakCenterMassHash' => \%peakCenterMassHash,
						'peakCenterMassIndexHash' => \%peakCenterMassIndexHash, 
						'peakIntensity' => \@peakIntensity, 
						'peakIndex' => \@peakIndex);
	return(\%msCentroidList, $intensityThreshold);
}

sub isMax {	## isMax determines whether the middle point reaches a local maximum
	my ($b2, $b1, $x, $a1, $a2) = @_;
	my $isMax = 0;
	if ($x > $b1 && $x > $a1) {
		$isMax = 1;		
	} elsif ($x == $b1 && $x > $b2 && $x > $a1) {
		$isMax = 1;
	} else {
		$isMax = 0;
	}
	return $isMax;
}

sub calcMinPeakIndex {	## calcMinPeakIndex determines the lower bound of a peak
	my ($ind, $intensity) = @_;
	while ($ind > 0 &&  $$intensity[$ind] != 0 && $$intensity[$ind - 1] <= $$intensity[$ind] ) {
		$ind--;
	}
	return $ind + 1;
}

sub calcMaxPeakIndex {	## calcMaxPeakIndex determines the upper bound of a peak
	my ($ind, $intensity) = @_;
	my $count = @$intensity;
	while ($ind < $count && $$intensity[$ind] != 0 && $$intensity[$ind + 1] <= $$intensity[$ind] ) { 
		$ind++;
	}
	return $ind - 1;
}

sub calcCenterMass {	## calcCenterMass estimates the mass (m/z) value at the center of a peak through least square fitting
	my ($nPoints, $minInd, $centerInd, $maxInd, $mass, $intensity) = @_;
	my $peakCenterMass;
	my $peakIntensity = 0;
	for (my $j = $minInd; $j <= $maxInd; $j++) {
		#$peakIntensity += $$intensity[$j];	# peakIntesntiy = Sum of intensities (or maximum intensity)
		if ($$intensity[$j] >= $peakIntensity) {
			$peakIntensity = $$intensity[$j];
		}
	}
	## There's a maximum point, but other points are zeros
	if ($minInd == $maxInd) {
		$peakCenterMass = $$intensity[$maxInd];
		return ($peakIntensity, $peakCenterMass);
	}
	## Right angled triangle-shaped peak
	if ($minInd == $centerInd) {
		$peakCenterMass = estimate2($$mass[$centerInd], $$mass[$centerInd + 1], $$intensity[$centerInd], $$intensity[$centerInd + 1]);
		return ($peakIntensity, $peakCenterMass);
	}
	## Left angled triangle-shaped peak
	if ($maxInd == $centerInd) {
		$peakCenterMass = estimate2($$mass[$centerInd - 1], $$mass[$centerInd], $$intensity[$centerInd - 1], $$intensity[$centerInd]);
		return ($peakIntensity, $peakCenterMass);
	}
	## Typical bell(triangle)-shaped peak
	if ($nPoints <= 3) {
		$peakCenterMass = estimate3($$mass[$centerInd - 1], $$mass[$centerInd], $$mass[$centerInd + 1], $$intensity[$centerInd - 1], $$intensity[$centerInd], $$intensity[$centerInd + 1]);
		return ($peakIntensity, $peakCenterMass);
	}
	
	## Typical bell(triangle)-shaped peak, use more than three data points	
	if ($nPoints > 3) {
		my $nleft;
		my $nright;
		my $d;
			
		## Select data points to be used in the least square fitting
		if ($nPoints % 2 == 1) {	## Simple case when $npoints is an odd number
			$d = int($nPoints / 2);
			$nleft = max($centerInd - $d, $minInd);
			$nright = min($centerInd + $d, $maxInd);
		} else {	## If $npoints is an even number, it would be a bit complicated
			$d = $nPoints / 2 - 1;
			$nleft = max($centerInd - $d, $minInd);
			$nright = min($centerInd + $d, $maxInd);
			if ($nleft != $minInd && $nright != $maxInd) {
				if ($$intensity[$nleft - 1] > $$intensity[$nright + 1]) {
					$nleft--;
				} else {
					$nright++;
				}
			} elsif ($nleft != $minInd) {
				$nleft--;
			} elsif ($nright != $maxInd) {
				$nright++;
			}
		}
		
		my @x;
		my @y;
		for (my $i = 0; $i < $nright - $nleft + 1; $i++){
			$x[$i] = $$mass[$nleft + $i];
			$y[$i] = $$intensity[$nleft + $i];
		}
		$peakCenterMass = estimateN(\@x, \@y, $centerInd - $nleft);
		return ($peakIntensity, $peakCenterMass);	
	}	
}

sub estimate2 {	## estimate2 gives a peakCenterMass by calculating an intensity-weighted average of mass
	my ($m1, $m2, $i1, $i2) = @_;
	my $peakCenterMass = ($m1 * $i1 + $m2 * $i2) / ($i1 + $i2);
	return $peakCenterMass;
}

sub estimate3 {	## estimate3 gives a peakCenterMass using the closed form of least-square fitting to a Gaussian-shaped peak
	my ($m1, $m2, $m3, $i1, $i2, $i3) = @_;
	my $l1 = log($i1);
	my $l2 = log($i2);
	my $l3 = log($i3);
	my $peakCenterMass = 0.5 * (($l1 - $l2) * ($m3 * $m3 - $m1 * $m1) - ($l1 - $l3) * ($m2 * $m2 - $m1 * $m1)) / (($l1 - $l2) * ($m3 - $m1) - ($l1 - $l3) * ($m2 - $m1));
	return $peakCenterMass;
}

sub estimateN {	## estimateN gives a peakCenterMass using least-square fitting to a Gaussian-shaped peak
	my ($x, $y, $ic) = @_;
	my @x = @{$x};
	my @y = @{$y};
	my $dm = $x[$ic];	## Intensity of centerInd
	my $regInput;
	for (my $i = 0; $i < @x; $i++) {
		$regInput -> [$i] = [log($y[$i]), ($x[$i]-$dm), ($x[$i]-$dm)**2];	## Mean-centering for x (singularity prevention???)		
	}
	my ($coeff, $Rsq) = linear_regression($regInput);	## Fitting to a Gaussian peak shape -> equivalent to the fitting of x and log(y) to a quadratic curve
	my $peakCenterMass = $dm - ($coeff -> [1]) / ($coeff -> [2]) * 0.5;
	return $peakCenterMass;
}

sub getMsHash {	
	my ($mzxml, $startScan, $endScan) = @_;	# Name of input .mzXML file
	my %ms;
	open (XML, "<", $mzxml) or die "Cannot open $mzxml file\n";
	my $indexOffset = getIndexOffset(*XML);
	my ($indexArray, $lastScan) = getIndexArray(*XML, $indexOffset);
	my $scanNumber = 0;
	my $scanIndex = 0;
	foreach (@$indexArray) {		
		$scanNumber++;
		if ($scanNumber < $startScan) {
			next;
		} elsif ($scanNumber > $endScan) {
			last;
		}
		print "  Gathering information of scan #$scanNumber out of $lastScan\r";
		my $index = $_;
		my $mslevel = getMsLevel(*XML, $index);
		next if ($mslevel != 1);
		my $rt = getRT(*XML, $index);
		my @peakArray;
		getPeak(*XML, \@peakArray, $index);	
		my $sizeArray = scalar(@peakArray);
		my @mz;
		my @int;
		for (my $i = 0; $i < $sizeArray; $i++) {
			$i++;
			push (@mz, $peakArray[$i - 1]);
			push (@int, $peakArray[$i]);
		}
		$ms{$scanIndex}{'RT'} = $rt;
		$ms{$scanIndex}{'scanNumber'} = $scanNumber;
		$ms{$scanIndex}{'mass'} = \@mz;
		$ms{$scanIndex}{'intensity'} = \@int;
		$scanIndex++;
	}
	close XML;
	return %ms;
}

sub getIndexOffset{
	(*XML) = @_;
	seek (XML, -120, 2);
	my ($indexOffset);
	while (<XML>) {
		next if (!/<indexOffset>/);
		chomp;
		$indexOffset = $_;
		$indexOffset =~ s/\s+<indexOffset>(\d+)<\/indexOffset>.*/$1/o;
		last;
	}
	return $indexOffset;
}

sub getIndexArray {	
	(*XML, my $indexOffset) = @_;
	my @indexArray;
	my $lastScan = 0;
	seek (XML, $indexOffset, 0);
	while (<XML>) {
		next if (/^\s+<index/);
		last if (/^\s+<\/index>.*/);
		chomp;
		next if (/scan/);
		$lastScan++;
		my $index = $_;
		$index =~ s/[\s\t]+\<offset id="\d+"[\s]*>(\d+)<\/offset>.*/$1/o;
		push (@indexArray, $index);
	}
	return (\@indexArray, $lastScan);
}

sub getRT {
	(*XML, my $scanIndex) = @_;
	my $rt;
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		next if (!/retentionTime=/);
		chomp;
		$rt = $_;
		$rt =~ s/.*retentionTime="PT([\d+\.]+)S".*/$1/o;
		last;
	}
	return $rt;
}

sub getScanNum {
	(*XML, my $scanIndex) = @_;
	seek (XML, $scanIndex, 0);
	my $scanNum;
	while (<XML>){
		next if (!/<scan\snum=\"\d+\"/);
		$_ =~ s/^M//g;
		chomp;
		$scanNum = $_;
		$scanNum =~ s/<scan\snum=\"(\d+)\"//;
		$scanNum = $1;
		last;
	}
	return $scanNum;
}

sub getMsLevel {
	(*XML, my $scanIndex) = @_;
	my $msLevel;	
	seek (XML, $scanIndex, 0);
	while (<XML>){
		next if (!/msLevel=/);
		chomp;
		$msLevel = $_;
		$msLevel =~ s/\s+msLevel="(\d)".*/$1/o;
		last;
	}
	return $msLevel;
}

sub getPeak{
	(*XML, my $peakArray, my $scanIndex) = @_;
	my ($peak, $peakLine);
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		if (/<peaks\sprecision="32"/ || /="m\/z-int"[\s]*>/){
			chomp;
			next if ($_ =~ /<peaks precision="32"\Z/);
			$peakLine = $_;
			if (/<peaks precision="32">/) {
				$peakLine =~ s/\s+<peaks precision="32"[\s\w\W\d\=\"]+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			} elsif (/contentType=\"m\/z\-int\"\>/) {
				$peakLine =~ s/\s+contentType=\"m\/z\-int\"\>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			} else {
				$peakLine =~ s/\s+="m\/z-int"[\s]*>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			}		
			last;
		} elsif (/\s+compressedLen="[\d]"\s+>([A-Za-z0-9\/\+\=]+)<\/peaks>/) {
			chomp;
			$peakLine = $_;
			$peakLine =~ s/\s+compressedLen="[\d]"\s+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			last;
		} else {
			next;
		}
	}
	#Base64 Decode
	$peak = decode_base64($peakLine);
	my @hostOrder32 = unpack("N*", $peak);
	for (@hostOrder32){
		my $float = unpack("f", pack("I", $_));
		push (@$peakArray, $float);
	}
}

sub indMax{
	my (@array) = @_;
	my $valMax = max(@array);
	my @indMax;
	for (my $i = 0; $i < scalar(@array); $i++) {
		if ($array[$i] == $valMax) {
			push @indMax, $i;
		}
	}
	return (@indMax);
}

sub makeTemporaryFeatureTable {
	my ($inputFile, $featureRef, $noisyScan, $subDir, $params) = @_;
	my @features = @{$featureRef};
#	my $tmpFile = basename($inputFile);
#	$tmpFile =~ s/\.mzXML//;
#	$tmpFile = $subDir . "/" . $tmpFile . ".tmp.feature";
	my $tmpFile = $subDir . ".tmp.feature";
	open (TMP, ">", $tmpFile);
	print TMP "index\tm\/z\tMS1 scan#\tIntensity\tS/N\tAssociated scans\tmin RT\tmax RT\tpercentage_true_peaks\n";		
	my ($minRT, $maxRT) = (0, 0);
	for (my $i = 0; $i < scalar(@features); $i++) {
		if ($minRT > min(@{$features[$i]{'scanRT'}})) {
			$minRT = min(@{$features[$i]{'scanRT'}});
		}
		if ($maxRT < max(@{$features[$i]{'scanRT'}})) {
			$maxRT = max(@{$features[$i]{'scanRT'}});
		}
	}
	my $index = 0;
	for (my $i = 0; $i < scalar(@features); $i++) {
		## Calculation of m/z for each feature
		my $nPoints = @{$features[$i]{'centerMz'}};
		my $numerator = 0;
		my $denominator = 0;
		for (my $j = 0; $j < $nPoints; $j++) {
			$numerator += $features[$i]{'centerMz'}[$j] * $features[$i]{'intensity'}[$j];
			$denominator += $features[$i]{'intensity'}[$j];
		}
		next if ($denominator == 0);	
		my $mz = $numerator / $denominator;	## Mass estimate calculated by the weighted mean of smoothed 3D-peak intensities
		$features[$i]{'3Dmass'} = $mz;
		
		## Filter out features with low intensity level
		my $strongIntensity = 0;
		my $strongScanNumber = 0;
		my %featureIntensity;	
		foreach (my $j = 0; $j < scalar(@{$features[$i]{'intensity'}}); $j++) {
			if ($strongIntensity < $features[$i]{'intensity'}[$j]) {
				$strongIntensity = $features[$i]{'intensity'}[$j];
				$strongScanNumber = $features[$i]{'scanNumber'}[$j];
			}
			$featureIntensity{$features[$i]{'intensity'}[$j]} = $features[$i]{'scanNumber'}[$j];
		}
		$$noisyScan{$strongScanNumber} = 500 if (!defined($$noisyScan{$strongScanNumber}));
		my $featureIntensityThreshold = $$params{'signal_noise_ratio'} * $$noisyScan{$strongScanNumber};
		if ($strongIntensity >= $featureIntensityThreshold) {
			my $minFeatureRT = min(@{$features[$i]{'scanRT'}});
			my $maxFeatureRT = max(@{$features[$i]{'scanRT'}});
			if (($maxFeatureRT - $minFeatureRT) / ($maxRT - $minRT) * 100 <= $$params{'max_percentage_RT_range'}) {
				$index++;
				print TMP $index, "\t", $mz, "\t", $strongScanNumber, "\t", $strongIntensity, "\t",
									$strongIntensity / $$noisyScan{$strongScanNumber}, "\t", $strongScanNumber, "\t",
									$minFeatureRT, "\t", $maxFeatureRT, "\t", 
									($maxFeatureRT - $minFeatureRT) / ($maxRT - $minRT) * 100, "\n";
			}
		} else {
			next;
		}
	}
	close (TMP);
	return ($tmpFile);
}

sub makeFeatureTable {
	my ($mzXML, $tmpFeatureFile) = @_;		
	open (TMP, "<", $tmpFeatureFile);
	my %hash;
	my %intensityHash;
	<TMP>;
	while (<TMP>) {
		chomp $_;
		my @data = split(/\t/, $_);
		my ($mz, $ms1ScanNum, $intensity) = ($data[1], $data[2], $data[3]);
		$hash{$ms1ScanNum}{$intensity}{$mz} = $_;
		$intensityHash{$mz} = $intensity;
	}
	close(TMP);
	
	## Investigate features from one with the highest intensity
	## If there's any feature produced from an isotopic peak of another feature,
	## those features will be merged 
	my %charge;
	my %isotope;
	foreach my $scan (keys %hash) {
		foreach my $intensity (sort {$b <=> $a} keys %{$hash{$scan}}) {
			foreach my $mz (keys %{$hash{$scan}{$intensity}}) {	
				my $min = 0;
				($scan - 50) < 0 ? $min = 0 : $min = ($scan - 50);				
				for (my $i = $min; $i < ($scan + 50); $i++) {
					my $chargeHash = findCharge($mz, \%{$hash{$i}});
					foreach my $chargedMz (keys %$chargeHash) {
						$charge{$scan}{$mz} = $$chargeHash{$chargedMz};
						if ($intensityHash{$mz} > $intensityHash{$chargedMz}) {
							$isotope{$i}{$chargedMz} = $$chargeHash{$chargedMz};				
						} else {
							$isotope{$scan}{$mz} = $$chargeHash{$mz};
						}
					}
				}
			}
		}
	}
	
	open (XML, "<", $mzXML) || die "Cannot open .mzXML file";
	my $indexOffset = getIndexOffset(*XML);
	my ($indexArray, $lastScan) = getIndexArray(*XML, $indexOffset);
	my %scan2RT;
	my %RT2scan;
	foreach my $index (@$indexArray) {
		my ($scan) = getScanNum(*XML, $index);	
		my ($rt) = getRT(*XML, $index);
		$scan2RT{$scan} = $rt;
		$RT2scan{$rt} = $scan;
	}
	close (XML);
	
	my $featureFile = $tmpFeatureFile;
	$featureFile =~ s/\.tmp\.feature/\.feature/;
	open (FEATURE, ">", $featureFile);
	print FEATURE "index", "\t", "m\/z", "\t", "z", "\t", "MS1scan#", "\t", "minMS1Scan#", "\t", "maxMS1Scan#", "\t", 
					"RT", "\t", "minRT", "\t", "maxRT", "\t", "Intensity", "\t", "S\/N", "\t", "Percentage of TF", "\n";
	my $index = 0;
	foreach my $scan (sort {$a <=> $b} keys %hash) {
		foreach my $intensity (sort {$b <=> $a} keys %{$hash{$scan}}) {
			foreach my $mz (sort {$a <=> $b} keys %{$hash{$scan}{$intensity}}) {	
				next if (defined($isotope{$scan}{$mz}));
				my @data = split("\t", $hash{$scan}{$intensity}{$mz});
				my ($SN, $minRT, $maxRT, $pctTruePeaks) = ($data[4], $data[6], $data[7], $data[8]);
				$index++;
				if (defined($charge{$scan}{$mz})) {
					print FEATURE $index, "\t", sprintf("%.12f", $mz), "\t", $charge{$scan}{$mz}, "\t", $scan, "\t", $RT2scan{$minRT}, "\t", $RT2scan{$maxRT}, "\t",
									sprintf("%.2f", $scan2RT{$scan}), "\t", $minRT, "\t", $maxRT, "\t",
									sprintf("%.0f", $intensity), "\t", sprintf("%.1f", $SN), "\t", sprintf("%.4f", $pctTruePeaks), "\n";
				} else {
					print FEATURE $index, "\t", sprintf("%.12f", $mz), "\t", "0", "\t", $scan, "\t", $RT2scan{$minRT}, "\t", $RT2scan{$maxRT}, "\t",
									sprintf("%.2f", $scan2RT{$scan}), "\t", $minRT, "\t", $maxRT, "\t",
									sprintf("%.0f", $intensity), "\t", sprintf("%.1f", $SN), "\t", sprintf("%.4f", $pctTruePeaks), "\n";
				}
			}
		}
	}
	close (FEATURE);
	return ($featureFile);
}

sub findCharge {
	my ($selectMz, $hash) = @_;
	my $C = 1.00335;
	my $intraPpm = 10;	## Decharge ppm
	my $maxCharge = 6;
	my ($lL, $uL) =  ($selectMz, $selectMz + $C + $selectMz * $intraPpm / 1e6); 
	my %chargeHash;	
	foreach my $intensity (sort {$b <=> $a} keys %$hash) {
		foreach my $mz (keys %{$$hash{$intensity}}) {
		# search the previous peak (only one peak) 
			if ($mz > $lL && $mz < $uL) {
				my $diff = 1 / abs($mz - $selectMz);
				my $roundDiff = sprintf("%.0f", $diff);
				next if ($roundDiff == 0 || $roundDiff > $maxCharge);				
				my $var = abs(abs($mz - $selectMz) - ($C / $roundDiff));
				next if ($var > ($intraPpm / 1e6) * $selectMz);   				
				my $charge = $roundDiff;
				$chargeHash{$mz} = $charge;
				return \%chargeHash;
			}
		}
	}
	return \%chargeHash;
}

sub median {
        my (@data) = @_;
        return unless @data;
        return $data[0] unless scalar(@data) > 1;
        @data= sort {$a <=> $b} @data;
        return $data[$#data / 2] if @data&1;
        my $mid= scalar(@data) /2;
        return ($data[$mid - 1] + $data[$mid]) / 2;
}
