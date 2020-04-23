#!/usr/bin/perl  

# Load all modules required by JUMPm 
## FindBin for getting current working directory
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";
use Spiders::Params;
use Spiders::ProcessingRAW;
use Spiders::ProcessingMzXML;
use Spiders::SGEQueue;
use Spiders::LSFQueue;
use Spiders::Path;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Basename;

####################
## Initialization ##
####################
my ($paramFile, $ms2Path) = @ARGV;
## Loading parameters
my $p = Spiders::Params -> new('-path' => $paramFile);
my $params = $p -> parseParams();
## Set the job management
my $queue;
if ($$params{'cluster'} == 1 && ($$params{'cluster_type'} eq "LSF" || $$params{'job_management_system'} eq "LSF")) {
	$queue = Spiders::LSFQueue->new();
} elsif ($$params{'cluster'} == 1 && ($$params{'cluster_type'} eq "SGE" || $$params{'job_management_system'} eq "SGE")) {
    $queue = Spiders::SGEQueue->new();
}

#########################################################################
## Database search to identify metabolites using individual .MS2 files ##
#########################################################################
my @ms2FileArray = glob("$ms2Path/*.MS2");
my $nFiles = scalar(@ms2FileArray);
my $maxJobs = 200;
my $filesPerJob;
if ($nFiles <= $maxJobs) {
	$maxJobs = $nFiles;
	$filesPerJob = 1;
} else {
	$filesPerJob = int($nFiles / $maxJobs) + 1;
}
my $nTotalJobs = int($nFiles / $filesPerJob - 0.0001) + 1;
my $nJobs = 0;
my $randNum = int(rand(100));
my %jobIDs = {};
for (my $i = 0; $i < $nTotalJobs; $i++) {
	$nJobs++;
	my $jobName = "sch_m_$nJobs";
	my $command = "";
	for (my $j = 0; $j < $filesPerJob; $j++) {
		my $k = $filesPerJob * $i + $j;
		last if ($k >= $nFiles);
		$command .= "perl $Bin/databaseSearch.pl $paramFile $ms2FileArray[$k]\n";
	}
	my $job = $queue -> submit_job($ms2Path, $jobName, $command);
	$jobIDs{$job} = 1;
	print "\r  $nJobs database search jobs are submitted";
}
print "\n  You submitted $nJobs jobs for database search\n";
checkJobStatus($nJobs, \%jobIDs, $queue);
print "\n";

#################
## Subroutines ##
#################
sub checkJobStatus {
	my ($nJobs, $jobIDs, $queue) = @_;
	my $nFinishedJobs = 0;
	my $jobInfo = 1;
	my ($username) = getpwuid($<);
	$| = 1;
	while($jobInfo) {
		my @jobStatusArray = @{$queue->get_running_jobs($jobIDs)};
		if (@jobStatusArray) {  # i.e. There are some jobs in the queue (may include other jobs like searching)
			my @jobIDsInQueue;
			foreach (@jobStatusArray) {
				my ($jobID) = ($_ =~ /([0-9]+)\s*/);
				if (defined $$jobIDs{$jobID}) {
					push (@jobIDsInQueue, $jobID);
				}
			}
			if (@jobIDsInQueue) { # i.e. There are jobs of interest in the queue
				$nFinishedJobs = $nJobs - scalar(@jobIDsInQueue);
				print "\r  $nFinishedJobs jobs are finished";
				if ($nFinishedJobs == $nJobs) {
					$jobInfo = 0;
				} else {
					$jobInfo = 1;
				}
			} else {        # i.e. There are jobs in the queue, but all jobs of interest are finished
				print "\r  $nJobs jobs are finished";
				$jobInfo = 0;
			}
		} else {        # i.e. There's no job in the queue
			print "\r  $nJobs jobs are finished";
			$jobInfo = 0;
		}
		sleep(5);
	}
	$| = 0;	
}
