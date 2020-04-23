#!/usr/bin/perl

package Spiders::Params;
use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = qw(error path parseParams  get_adducts parameterCheck get_full_text  write_text  reconstruct  rewrite  update  get  set_default_dir  get_defaults); 
my $ERROR;

sub new {
	my ($class, %args) = @_;
	my $self = {};
	$self->{PATH} = $args{-path};
	bless $self;
#	$self->parseParams();
#	$ERROR = $@ if($@);
	return $self;
}

sub error {
	return $ERROR;
}

sub path {
	my $self = shift;
	$self -> {PATH} = shift if (@_);
	return $self -> {PATH};
}
	
sub parseParams {
	my ($self) = @_;  
	my ($paramhash, $paramfile) = @_;	
	my %parameters;
	## Default parameter settings
	$parameters{'output_name'} = 1;
	$parameters{'labeled_data'} = 1;
	$parameters{'labeled_ID_method'} = 1;
	$parameters{'mode'} = 1;	
	$parameters{'unlabel_scan'} = 1;
	$parameters{'unlabel_min_mass'} = 1;
	$parameters{'unlabel_max_mass'} = 1;
	$parameters{'database'} = 1;
	$parameters{'mass_formula_database'} = 1;
	$parameters{'structure_database'} = 1;
	$parameters{'formula_mass_tolerance_searching'} = 1;
	$parameters{'formula_mass_tolerance_pairing'} = 1;
	$parameters{'c12_c13_tolerance'} = 1;
	$parameters{'c12_n15_tolerance'} = 1;
	$parameters{'relative_isotopes_intensity'} = 1;
	$parameters{'min_pair_correlation'} = 1;
	$parameters{'cluster_tolerance'} = 1;
	$parameters{'first_scan_extraction'} = 1;
	$parameters{'last_scan_extraction'} = 1;
	$parameters{'tolerance_unit'} = 1;	
	$parameters{'isolation_window'} = 1;
	$parameters{'mass_correction'} = 1;
	$parameters{'loading_normalization'} = 1;
	$parameters{'decharge_ppm'} = 1;
	$parameters{'deisotope_ppm'} = 1;

	$parameters{'percentage_ms2_peaks'} = 1;
	$parameters{'frag_mass_tolerance'} = 1;		
	$parameters{'frag_mass_tolerance_unit'} = 1;
	$parameters{'matched_scan_dist'} = 1;
	
	$parameters{'mass_shift'} = 1;
	$parameters{'data_acquisition_mode'} = 1;
	$parameters{'intensity_threshold'} = 1;
	$parameters{'intensityThresholdScale2D'} = 1;
	$parameters{'noise_detection_method'} = 1;
	$parameters{'noisePercentage2D'} = 1;
	$parameters{'noisePercentage3D'} = 1;
	$parameters{'processor_number'} = 1;
	$parameters{'fragment_method'} = 1;
	$parameters{'cluster'} = 1;
	$parameters{'job_management_system'} = 1;
	$parameters{'signal_noise_ratio'} = 1;
	$parameters{'adduct'} = 1;
	$parameters{'library_search'} = 1;
	$parameters{'decoy_strategy'} = 1;
	$parameters{'max_percentage_RT_range'} = 1;
	$parameters{'min_peak_intensity'} = 1;
	$parameters{'skipping_scans'} = 1;
	$parameters{'mass_tolerance_peak_matching'} = 1;   

	## Parameter file loading
	my $path = $self -> {PATH};
	if (open (P, "<", "$path")) {
		my $line;
		my $hash = {};    
		while ($line = <P>) {
			my $linehash = {};
			my $comments = "";      
			if ($line =~ s/\s*([;\#].*)$//) {
				$comments = $1;
			}      
			$linehash -> {Comments} = $comments;
			if ($line =~ /^(.+?)\s*=\s*(.+)$/) {
				my ($key, $data) = ($1, $2);
				$data =~ s/\s+$//o;
				if ($key eq "feature_files") {
					push (@{$hash -> {$key}}, $data);
				} else {
					$hash -> {$key} = $data;
				}
			} elsif ($line !~ /=/ && $line !~ /\#/ && $line =~ /\.feature/ && defined($hash -> {'feature_files'})) {
				chomp ($line);
				push (@{$hash -> {'feature_files'}}, $line);
			}
		}
		close P;
		$self -> {PARAMETERS} = $hash;
		return $hash;
	}
	return 0;
}

sub parameterCheck {
	my ($self, $params) = @_;
	if (($params -> {last_scan_extraction} - $params -> {first_scan_extraction}) < 50) {
		print "  Please set a larger scan region (e.g. >100) for peak dection\n";
		print "  For unlabeled data, the searched scan is defined by a parameter of unlabel_scan instead of these two parameters\n";
		print "  The scan region has to cover unlabel_scan\n\n";		
		exit;
	}
	if (!-e($params -> {mass_formula_database})) {
		print "  The specified mass formula database can not be found\n";
		exit;
	}
	if (!-e($params -> {structure_database})) {
		print "  The specified structure database can not be found\n";
		exit;
	}	
}

=head
sub get_adducts {
	my ($self, $parameter) = @_;
	my $adducts;
	foreach my $param (keys %$parameter) {
		if ($param =~ /adduct\_(\w+)/) {
			$adducts -> {$1} = $parameter -> {$param};
		}
	}
	return $adducts;
}


sub get_full_text {
	my ($self) = @_;
	local $/ = undef;
	my $path = $self -> {PATH} || return;
	open (P, "<", "$path") || return;
	my $p = <P>;
	close P;
	return $p;
}

sub write_text {
	my ($self) = @_;
}

sub reconstruct {
	my ($self, $newpath) = @_;
	my $path = $self -> {PATH};
	open (P, "<", "$path") || die "Couldn't open $path: $!";
	my @lines = <P>;
	close P;
  	
  	my $mask = "%-40s";
	open (NEW, ">",  "$newpath") || die "Couldnt't open $newpath: $!";
	my $line;
	my @newlines = ();
	my $parameters = $self -> {PARAMETERS};
	foreach $line (@lines){
		if ($line =~ /^(.+?)\s*=\s*(.+)$/ && exists $parameters -> {$1}) {
			printf (NEW $mask, "$1 = $parameters -> {$1} -> {data}");
			print NEW $parameters -> {$1} -> {Comments} . "\n";
		} else {
			print NEW $line;
		}
	}
	close NEW;
}

sub reconstruct_new {
	my ($self, $newpath) = @_;
	my $path = $self -> {PATH};
	open (P, "<", "$path") || die "Couldn't open $path: $!";
	my @lines = <P>;
	close P;

	my $mask = "%-40s";
	open (NEW, ">", "$newpath") || die "Couldnt't open $newpath: $!";
	my $line;
	my @newlines = ();
	my $parameters = $self -> {PARAMETERS};
	foreach $line (@lines){
		next if ($line =~ /duc/);
		if ($line =~ /SEQUEST_ENZYME_INFO/) {
			@lines = ();
		}
		if ($line =~ /^(.+?)\s*=\s*(.+)$/ && exists $parameters -> {$1}) {
			printf (NEW $mask, "$1 = $parameters -> {$1} -> {data}");
			print NEW $parameters -> {$1} -> {Comments}."\n";
		} else {
			print NEW $line;
		}
	}
	close NEW;
}

sub rewrite {
	my ($self) = @_;
	return $self -> reconstruct($self -> {PATH});
}

sub update {
	my ($self, $key, $newval) = @_;
	my $parameters = $self -> {PARAMETERS};
	if (exists $parameters -> {$key}){
		$parameters -> {$key} -> {data} = $newval;
	}
}

sub get {
	my ($self, $key) = @_;
	my $parameters = $self -> {PARAMETERS};
	if ($key eq 'current_enzyme') {
		return $self -> current_enzyme();
	}
	if (my $p = $parameters->{$key}) {
		return $p -> {data};
	}
	return 'Unknown';
}

sub set_default_dir {
	my ($self, $defaultdir) = @_;
	if (opendir (DIR, $defaultdir)){
		$defaultdir =~ s/\/$//o;
		$defaultdir .= '/';
		$self -> {defaultdir} = $defaultdir;
		my @defaults = map {"$defaultdir" . $_} grep {/\.params/} readdir(DIR);
		$self -> {defaults} = \@defaults;
	}
}

sub get_defaults {
	my ($self) = @_;
	return $self -> {defaults};
}
=cut

1;