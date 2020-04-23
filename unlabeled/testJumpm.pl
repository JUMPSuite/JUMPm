#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#my $nFiles = 4200;
#my $maxJobs = 200;
#my $filesPerJob;
#if ($nFiles <= $maxJobs) {
#	$maxJobs = $nFiles;
#	$filesPerJob = 1;
#} else {
#	$filesPerJob = int($nFiles / $maxJobs) + 1;
#}
#my $nTotalJobs = int($nFiles / $filesPerJob - 0.0001) + 1;
#my $nJobs = 0;
#for (my $i = 0; $i < $nTotalJobs; $i++) {
#	$nJobs++;		
#	for (my $j = 0; $j < $filesPerJob; $j++) {
#		my $k = $filesPerJob * $i + $j;
#		if ($k >= $nFiles) {
#			last;
#		}
#		print "job# $nJobs, file# $k\n";
#	}
#	print "$i\n";
#}

#my $paramFile;
#GetOptions('-p=s' => \$paramFile,);
#if (!-e ($paramFile)) {
#	print "Please input the parameter file\n\n";
#	exit;
#}
#print "$paramFile\n";
#print Dumper(@ARGV), "\n";

#my $str = "[M-3H]3-";
#my $charge = 3;
#my $coeff = $charge + 1;
#$str =~ s/M-\d/M-$coeff/;
#print "$str\n";

my $str = "[M-2H]2-";
my ($charge) = $str =~ /\[*.\](\d)[+-]/;
my $coeff = $charge + 1;
$str =~ s/(M-.*)H/M-\Q$coeff\EH/;
#$str =~ s/(M[+-].*H)/$1+Na/;
print "$str\n";

#my $str = "[M-2H]2-";
#my $adduct = "HCOO";
#$str =~ s/(M[+-].*H)/$1+$adduct/;
#print "$str\n";