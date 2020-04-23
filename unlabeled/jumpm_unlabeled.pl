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
use Parallel::ForkManager;

####################
## Initialization ##
####################
my $library = $Bin;
my $RLibrary = $Bin . "/R/";
if (!-e($RLibrary)) {
	print "R scripts are NOT accessible from jump -m\n";
	exit;
}
my $paramFile;
GetOptions('-p=s' => \$paramFile,);
if (!-e ($paramFile)) {
	print "Please input the parameter file\n\n";
	exit;
}
$paramFile = abs_path($paramFile);

## Define a log file usings time() function
my $tmpLog = time();
my $LOG;
## For unlabeled data, the log file should be open with append option to include R script captions
open ($LOG, ">>", "./$tmpLog");
print "\n\n  Initializing JUMPm program\n\n";
print $LOG "\n\n  Initializing JUMPm program\n\n";
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 

## Loading parameters
my $p = Spiders::Params -> new('-path' => $paramFile);
my $params = $p -> parseParams();

## Set the job management
my $queue;
my $MAX_PROCESSES = 4;
if ($$params{'cluster'} == 1 && $$params{'job_management_system'} eq "LSF") {
	$queue = Spiders::LSFQueue->new();
} elsif ($$params{'cluster'} == 1 && $$params{'job_management_system'} eq "SGE") {
	$queue = Spiders::SGEQueue->new();
} elsif ($$params{'cluster'} == 0) {
	$queue = new Parallel::ForkManager($MAX_PROCESSES);
}

##################
## Main routine ##
##################
my $subDir;
my %fileHash;
my @fileArray;
my @featureFileArray;
if ($$params{'skip_feature_detection'} == 0) {
	## Perform feature detection using input .raw (or .mzXML) files	
	## Preparation of directories according to the names of input files
	print "  Using the following raw files:\n";
	print $LOG "  Using the following raw files:\n";
	foreach my $arg (sort @ARGV) {
		if ($arg =~ /.[raw|RAW|mzXML]/) {
			print "  $arg\n";
			print $LOG "  $arg\n";
		}
	}
	my $nJobs = 0;
	my $randNum = int(rand(100));
	my %jobIDs;
	foreach my $arg (sort @ARGV) {
		#############################################################
		## Handling .raw (.mzXML) files and generating directories ##
		#############################################################
		$arg = basename($arg);	
		my @extensionList = (".raw", ".RAW", ".mzXML");
		if ($arg =~ /.[raw|RAW|mzXML]/) {
			push (@fileArray, $arg);
			my ($fileName, $directory, $extension) = fileparse($arg, @extensionList);
			## If $arg = (/Example/unlabeled_data/)12C_HILIC_Neg.mzXML, then
			## $filename = "12C_HILIC_Neg"
			## $directory = "/Example/unlabeled_data/"
			## $extension = "mzXML"
			my $dataDir = "$directory/$fileName";
			$dataDir =~ s/\/\//\//g;	## Change any double slashes (//) to single ones (/)
			system (qq(mkdir $dataDir >/dev/null 2>&1));
			system (qq(mv $arg $dataDir >/dev/null 2>&1));
			## At this point, $dataDir = /Example/unlabeled_data/12C_HILIC_Neg
			## Input file will be moved to the above $dataDir, /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.mzXML

			if (!-e($dataDir)) {
				print "  There is no such file ($arg)!\n";
				exit;
			}
			my $path = new Spiders::Path($dataDir);
			my $list = $path -> makeDirectoryList();
			if (@$list) {
				$subDir = $path -> chooseDirEntry($list, "  Choose a .out file directory", $subDir);
			} else {
				$subDir	= $fileName . ".1";
			}
			print "  Using: $subDir\n";
			print $LOG "  Using: $subDir\n";		
			$path -> addSubdir($subDir);
			$subDir =  $path -> baseDir() . "/$subDir";	## e.g. /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.1/		
			my $inputFile = "$dataDir/$arg";
			$fileHash{$inputFile} = $subDir;
			
			####################################################
			## File conversion (.raw to .mzXML), if necessary ##
			####################################################
			if ($inputFile =~ /\.mzXML/) {
				$inputFile =~ s/\.mzXML//g;
			}
			my $procRaw = new Spiders::ProcessingRAW();
			$procRaw -> setLogFile($LOG);		
			$procRaw -> setRawFile($inputFile);
			my $procXml = new Spiders::ProcessingMzXML();	
			$procXml -> setDtaPath($subDir);
			$procXml -> setLogFile($LOG);		
			print "  Handling input file(s)\n";
			print $LOG "  Handling input file(s)\n";
			my $mzXML;
			if ($$params{'data_acquisition_mode'} == 1) {	## Centroid mode
				$mzXML = $procRaw -> raw2mzXML("centroid");
			} elsif ($$params{'data_acquisition_mode'} == 2) {
				$mzXML = $procRaw -> raw2mzXML("profile");
			} else {
				system (qq(cp $tmpLog $subDir/JUMPm.log));		
				print "Please set the proper 'data_acquisition_mode' parameter\n";
				exit;			
			}
	
			#######################
			## Feature detection ##
			#######################
			$nJobs++;
			$mzXML = abs_path($mzXML);
			$dataDir = abs_path($dataDir);
			my $jobName = "feat_m_$nJobs";
			my $command = "perl $Bin/featureDetection.pl $paramFile $mzXML $subDir";
			my $featureFile = $subDir . "\.feature";
			push (@featureFileArray, $featureFile);
			if ($$params{'cluster'} == 0) {
				$queue -> start and next;
				my $file = File::Spec->catfile($dataDir, "$jobName.sh");
				open (JOB, ">", $file) or die "Cannot create a job file\n";
				print JOB "#!/bin/bash\n";
				print JOB "#\$ -S /bin/bash\n";
				print JOB "#\$ -N $jobName\n";
				print JOB "#\$ -e $dataDir/$jobName.e\n";
				print JOB "#\$ -o $dataDir/$jobName.o\n";
				print JOB $command;
				close (JOB);
				system ("cd $dataDir && sh '$jobName.sh' >/dev/null 2>&1");
				$queue -> finish; # Terminates the child process
			} elsif ($$params{'cluster'} == 1) {
				my $job = $queue -> submit_job($dataDir, $jobName, $command);
				$jobIDs{$job} = 1;
				print "\r    $nJobs feature detection job(s) is/are submitted";
			}
		}
	}
	if ($$params{'cluster'} == 0) {
		print "\n  You submitted $nJobs job(s) for feature detection\n";
		print "\n  It may take a while to finish those jobs in a single server. Please be patient...\n";
		$queue -> wait_all_children;
	} elsif ($$params{'cluster'} == 1) {
		print "\n    You submitted $nJobs job(s) for feature detection\n";
		checkJobStatus($nJobs, \%jobIDs, $queue);
	}
} elsif ($$params{'skip_feature_detection'} == 1) {
	print "  Feature detection step is skipped\n";
	print "  The following feature file(s) is/are used\n";
	print $LOG "  Feature detection step is skipped\n";
	print $LOG "  The following feature file(s) is/are used\n";
	## Skip feature detection
	@featureFileArray = @{$$params{'feature_files'}};
	foreach my $featureFile (sort {$a cmp $b} @featureFileArray) {
		print "  $featureFile\n";
		print $LOG "  $featureFile\n";
		my $key = $featureFile;
		my $value = $featureFile;
		$key =~ s/\d+\.feature/mzXML/;
		$value =~ s/\.feature//;
		$fileHash{$key} = $value;
		push (@fileArray, basename($key));
	}
	$subDir = $featureFileArray[0];
	$subDir =~ s/\.feature$//;
} else {
	print "  'skip_feature_detection' parameter should be set properly\n";
	print "  Please check the parameter and try again\n";
	exit;
}
## At this point, feature files are generated and located in the subdirectories, e.g.
## current working directory = /Example/unlabeled_data/
## /Example/unlabeled_data/jump_m.params
## /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.mzXML
## /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.mzXML
## /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.1.tmp.feature
## /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.1.feature
## /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.1/ ...
## ...
## $dataDir = /Example/unlabeled_data/12C_HILIC_Neg
## $subDir =  /Example/unlabeled_data/12C_HILIC_Neg/12C_HILIC_Neg.1/
print "\n";
print $LOG "\n";

#######################
## Feature alignment ##
#######################
=head
## Create alignment.params file
open (AP, ">", "alignment.params");
foreach my $key (sort {$a cmp $b} keys %$params) {
	if ($key eq "tol_initial") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "sd_width") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "rescue") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "pct_full_alignment") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "output_name") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "rt_tolerance_unit") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "rt_tolerance_value") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "mz_tolerance_unit") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "mz_tolerance_value") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "reference_feature") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "skip_loading_bias_correction") {
		print AP "$key" . " = " . $$params{$key} . "\n";
	}
}

## Execute R script for feature alignment
my $featureFiles = join(",", @featureFileArray);
my $alignDir = "./align_" . $$params{'output_name'};
system(qq(mkdir $alignDir >/dev/null 2>&1));
my $command = "Rscript $Bin/R/alignment.R alignment.params $featureFiles $tmpLog $Bin $alignDir";
system ($command);
## $alignDir = ./align_test/ (when $$params{'output_name'} = test)
## abs_path($alignDir) = /Example/unlabeled_data/align_test/
=cut

my $featureFiles = join(",", @featureFileArray);
my $alignDir = "./align_" . $$params{'output_name'};
system(qq(mkdir $alignDir >/dev/null 2>&1));
my $command = "Rscript $Bin/R/alignment.R $paramFile $featureFiles $tmpLog $Bin $alignDir";
system ($command);

#########################################
## Processing MS2 spectra for features ##
#########################################
=head
## Create featureToMs2.params file
open (FP, ">", "featureToMs2.params");
my $fullyAlignedFeatureFile = abs_path($alignDir) . "/" . $$params{'output_name'} . "_fully_aligned.feature";
print FP "feature_file = $fullyAlignedFeatureFile\n";
## "tol_isolation" in featureToMs2 should be the same as "isolation_window" parameter
my $tol_isolation = $$params{'isolation_window'} / 2;
print FP "tol_isolation = $tol_isolation\n";
foreach my $key (sort {$a cmp $b} keys %$params) {
	if ($key eq "tol_precursor") {
		print FP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "tol_intra_ms2_consolidation") {
		print FP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "tol_inter_ms2_consolidation") {
		print FP "$key" . " = " . $$params{$key} . "\n";
	} elsif ($key eq "mode") {
		print FP "$key" . " = " . $$params{$key} . "\n";
	}
}
for (my $i = 0; $i < scalar(@fileArray); $i++) {
	foreach my $key (keys %fileHash) {
		if ($fileArray[$i] eq basename($key)) {
			my $mzXML = $key;
			$mzXML =~ s/\.(raw|RAW|Raw)/\.mzXML/;
			if ($i == 0) {
				print FP "mzXML = " . $mzXML . "\n";
			} else {
				print FP "$mzXML\n";
			}
		}
	}
}
close (FP);

## Execute R script for processing MS2 spectra
my ($ms2Path) = $subDir =~ /(\.\d+$)/;
$ms2Path = abs_path($alignDir) . "/align_" . $$params{'output_name'} . $ms2Path;
$command = "Rscript $Bin/R/featureToMs2.R featureToMs2.params $ms2Path $tmpLog";
system ($command);
## $ms2Path = /Example/unlabeled_data/align_test/align_test.1/
## i.e. $ms2Path indicates the path where .MS2 files for fully-aligned features are located
print "\n";
print $LOG "\n";
=cut

my $fullyAlignedFeatureFile = abs_path($alignDir) . "/." . $$params{'output_name'} . "_fully_aligned.feature";	## Load "hidden" fully aligned feature file
my @tmpArray;
for (my $i = 0; $i < scalar(@fileArray); $i++) {
	foreach my $key (keys %fileHash) {
		if ($fileArray[$i] eq basename($key)) {
			my $mzXML = $key;
			$mzXML =~ s/\.(raw|RAW|Raw)/\.mzXML/;
			push (@tmpArray, $mzXML);
		}
	}
}
my $mzXMLs = join(",", @tmpArray);
my ($ms2Path) = $subDir =~ /(\.\d+$)/;
$ms2Path = abs_path($alignDir) . "/align_" . $$params{'output_name'} . $ms2Path;
$command = "Rscript $Bin/R/featureToMs2.R $paramFile $fullyAlignedFeatureFile $mzXMLs $ms2Path $tmpLog";
system ($command);

##================##
## Library search ##
##================##
if ($$params{'library_search'} == 1) {	
	$command = "perl $Bin/librarySearch.pl $paramFile $fullyAlignedFeatureFile $ms2Path\n";
	system ($command);
} else {
	print "  Libraray search is skipped according to the parameter\n";
	print $LOG "  Libraray search is skipped according to the parameter\n";
}
print "\n";
print $LOG "\n";

##=================##
## Database search ##
##=================##
if ($$params{'database_search'} == 1) {
	#########################################################################
	## Database search to identify metabolites using individual .MS2 files ##
	#########################################################################
	print "  Database search for the identification of metabolites\n";
	print "  =====================================================\n\n";
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
	my %jobIDs;
	if ($$params{'cluster'} == 0) {
		$queue = new Parallel::ForkManager($MAX_PROCESSES);
		print "  Many jobs will be processed in a single server\n";
		print "  It may take a while. Please be patient ... \n";
	}
	for (my $i = 0; $i < $nTotalJobs; $i++) {
		$nJobs++;
		my $jobName = "sch_m_$nJobs";
		my $command = "";
		for (my $j = 0; $j < $filesPerJob; $j++) {
			my $k = $filesPerJob * $i + $j;
			last if ($k >= $nFiles);
			$command .= "perl $Bin/databaseSearch.pl $paramFile $ms2FileArray[$k]\n";
		}
		if ($$params{'cluster'} == 0) {
			$queue -> start and next;
			my $file = File::Spec->catfile($ms2Path, "$jobName.sh");
			open (JOB, ">", $file) or die "Cannot create a job file\n";
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -S /bin/bash\n";
			print JOB "#\$ -N $jobName\n";
			print JOB "#\$ -e $ms2Path/$jobName.e\n";
			print JOB "#\$ -o $ms2Path/$jobName.o\n";
			print JOB $command;
			close (JOB);
			system ("cd $ms2Path && sh '$jobName.sh' >/dev/null 2>&1");
			$queue -> finish; # Terminates the child process
		} elsif ($$params{'cluster'} == 1) {
			my $job = $queue -> submit_job($ms2Path, $jobName, $command);
			$jobIDs{$job} = 1;
			print "\r    $nJobs database search jobs are submitted";
		}
	}
	if ($$params{'cluster'} == 0) {
		$queue -> wait_all_children;
	} elsif ($$params{'cluster'} == 1) {
		print "\n    You submitted $nJobs jobs for database search\n";
		checkJobStatus($nJobs, \%jobIDs, $queue);
	}
	print "\n\n";
	
	######################################################
	## Calculation of scores for identified metabolites ##
	######################################################
	my @ms2OutFileArray = glob("$ms2Path/*.MS2.out");
	$nFiles = scalar(@ms2OutFileArray);
	if ($nFiles <= $maxJobs) {
		$maxJobs = $nFiles;
		$filesPerJob = 1;
	} else {
		$filesPerJob = int($nFiles / $maxJobs) + 1;
	}
	$nTotalJobs = int($nFiles / $filesPerJob - 0.0001) + 1;
	$nJobs = 0;
	$randNum = int(rand(100));
	%jobIDs = ();
	if ($$params{'cluster'} == 0) {
		$queue = new Parallel::ForkManager($MAX_PROCESSES);
		print "  Many jobs to calculate scores of identified metabolites will be processed in a single server\n";
		print "  It may take a while. Please be patient ... \n";
	}
	for (my $i = 0; $i < $nTotalJobs; $i++) {
		$nJobs++;
		my $jobName = "score_m_$nJobs";
		my $command = "";
		for (my $j = 0; $j < $filesPerJob; $j++) {
			my $k = $filesPerJob * $i + $j;
			last if ($k >= $nFiles);
			$command .= "perl $Bin/calculateScores.pl $paramFile $ms2OutFileArray[$k]\n";
		}
                if ($$params{'cluster'} == 0) {
                        $queue -> start and next;
                        my $file = File::Spec->catfile($ms2Path, "$jobName.sh");
                        open (JOB, ">", $file) or die "Cannot create a job file\n";
                        print JOB "#!/bin/bash\n";
                        print JOB "#\$ -S /bin/bash\n";
                        print JOB "#\$ -N $jobName\n";
                        print JOB "#\$ -e $ms2Path/$jobName.e\n";
                        print JOB "#\$ -o $ms2Path/$jobName.o\n";
                        print JOB $command;
                        close (JOB);
                        system ("cd $ms2Path && sh '$jobName.sh' >/dev/null 2>&1");
                        $queue -> finish; # Terminates the child process
		} elsif ($$params{'cluster'} == 1) {
			my $job = $queue -> submit_job($ms2Path, $jobName, $command);
			$jobIDs{$job} = 1;
			print "\r    $nJobs score calculation jobs are submitted";
		}
	}
	if ($$params{'cluster'} == 0) {
		$queue -> wait_all_children;
	} elsif ($$params{'cluster'} == 1) {
		print "\n    You submitted $nJobs jobs for database search\n";
		checkJobStatus($nJobs, \%jobIDs, $queue);
	}
	print "\n";
	
	##################################
	## Create a result table (file) ##
	##################################
	my $fullyAlignedFeatureFile = abs_path($alignDir) . "/" . $$params{'output_name'} . "_fully_aligned.feature";	## Load "new" fully aligned feature file
	my $command = "perl $Bin/generateResult.pl $ms2Path $fullyAlignedFeatureFile $paramFile";
	system ($command);

=head
	print "\n  Tables containing results are being generated\n";
	$ms2Path =~ s/\/$//;
	my $outFile1 = $ms2Path . ".spectrum_matches";  ## .spectrum_matches file containing all target and decoy structures with mscores for each feature
	my $outFile2 = $ms2Path . ".structure"; ## .structure file containing a unique (target or decoy structure
	open (OUT1, ">", $outFile1) or die "Cannot open $outFile1\n";
	open (OUT2, ">", $outFile2) or die "Cannot open $outFile2\n";
	my $header = "FeatureNo\tMass\tFormula\tName\tStructure\tInChi\tType\tAdduct\t" .
	                        "NumMeasuredPeaks\tNumTheoreticalPeaks\tNumMatchedPeaks\tMscore\tIntensityRatio\n";
	print OUT1 $header;
	print OUT2 $header;
	
	## Sort .score files according to the feature number
	my @scoreFiles = glob("$ms2Path/*.score");
	my @featureNums;
	foreach my $scoreFile (@scoreFiles) {
		my ($featureNum) = $scoreFile =~ /f(\d+)\.MS2\.score/;
		push (@featureNums, $featureNum);
	}
	my @index = sort {$featureNums[$a] <=> $featureNums[$b]} (0..$#featureNums);
	@featureNums = @featureNums[@index];
	@scoreFiles = @scoreFiles[@index];
	
	## Read each .score file and write to output files
	for (my $i = 0; $i < scalar(@scoreFiles); $i++) {
		open (SCORE, "<", $scoreFiles[$i]) or die "Cannot open $scoreFiles[$i]\n";
		my $bestEntry;  ## The best structure which has the highest "mscore"
		my $bestMscore = 0;
		while (<SCORE>) {
			chomp ($_);
			my $line = $_;
			next if ($line =~ "Index");
	
			## Replace the first element of $line (i.e. contents of .score file) with 'feature number' 
			my @elems = split(/\t/, $line);
			$elems[0] = $featureNums[$i];
			$line = join("\t", @elems);
			
			## Write to .spectrum_matches file
			print OUT1 $line, "\n";
	
			## Choose the entry with the highest mscore (i.e. best entry)
			my $mscore = (split(/\t/, $line))[-2];  ## mscore is the second last element
			if ($mscore == -0) {
				$mscore = 0;
			}
			if (!defined ($bestEntry)) {
				$bestEntry = $line;
				$bestMscore = $mscore;
			} else {
				if ($mscore > $bestMscore) {
					$bestEntry = $line;
					$bestMscore = $mscore;
				}
			}
		}
		close (SCORE);
		print OUT2 $bestEntry, "\n";
	}
	close (OUT1);
	close (OUT2);	
=cut
} else {
	print "  Database search is skipped according to the parameter\n";
	print $LOG "  Database search is skipped according to the parameter\n";	
}


#######################################
## Copy parameter files and LOG file ##
#######################################
system ("cp $paramFile $alignDir");
#system ("mv alignment.params $alignDir");
#system ("mv featureToMs2.params $alignDir");
system ("cp $tmpLog $alignDir/JUMPm.log");
system ("rm $tmpLog");

print "\n  Jump -m for unlabeled data is finished\n\n";

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
				print "\r    $nFinishedJobs jobs are finished";
				if ($nFinishedJobs == $nJobs) {
					$jobInfo = 0;
				} else {
					$jobInfo = 1;
				}
			} else {        # i.e. There are jobs in the queue, but all jobs of interest are finished
				print "\r    $nJobs jobs are finished";
				$jobInfo = 0;
			}
		} else {        # i.e. There's no job in the queue
			print "\r    $nJobs jobs are finished";
			$jobInfo = 0;
		}
		sleep(5);
	}
	$| = 0;	
}
