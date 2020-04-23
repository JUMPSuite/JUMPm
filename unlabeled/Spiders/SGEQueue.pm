package Spiders::SGEQueue;
use strict;
use File::Spec;

sub new {
	my ($class, %arg) = @_;
	my $self = {};
	bless $self, $class;
	return $self;
}

sub submit_job {
	my ($self, $path, $name, $cmd) = @_;
	my $file = File::Spec->catfile($path, "$name.sh");
	open (JOB, ">", $file) or die "Cannot create a job file\n";
	print JOB "#!/bin/bash\n";
	print JOB "#\$ -S /bin/bash\n";
	print JOB "#\$ -N $name\n";
	print JOB "#\$ -e $path/$name.e\n";
	print JOB "#\$ -o $path/$name.o\n";
	print JOB $cmd;
	close (JOB);
#	my $command = qq(cd $path && qsub -cwd -pe mpi 8 -l mem_free=8G,h_vmem=8G $file);
	my $command = qq(cd $path && qsub -cwd -pe mpi 8 $file);
	my $job = lc(qx[$command]);
	chomp ($job);
	if ($job =~ /job (\d+)/) {
		$job = $1;
	} else {
		warn "could not parse qsub output";
	}
	return $job;
}

sub get_running_jobs {
	my ($self, $jobshash) = @_;
	my $command = "qstat";
	my $jobStatus = qx[$command];
	my @lines = split(/\n/, $jobStatus);
	shift @lines;
	my @retval;
	foreach my $l (@lines) {
		my @toks = split(/\s+/, $l);
		## Sometimes, $toks[0] is a blank. In this case, shift the memory @toks (2019/10/15)
		if ($toks[0] eq '') {
			shift @toks;
		}
		if (defined($jobshash->{$toks[0]})) {
			push (@retval, $toks[0]);
		}
	}
	return \@retval;
}
1;
