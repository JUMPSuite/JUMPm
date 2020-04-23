#!/usr/bin/perl

package Spiders::ProcessingRAW;
use strict;
use warnings;
use File::Basename;
use vars qw($VERSION @ISA @EXPORT);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = qw(set_raw_file get_raw_file get_rawfile_basename set_log_file get_log_file get_rawfile_dirname set_converter get_converter raw2mzXML get_mzXML_file);
 
sub new {
	my ($class, %arg)=@_;
    my $self = {
        _rawFile => undef,
		_mzXMLFile => undef,
		_converter => undef,
    };
    bless $self, $class;
	return $self;
}

sub setRawFile {
	my ($self, $rawFile) = @_;
	$self -> {_rawFile} = $rawFile;
}

sub setLogFile {
	my ($self, $log) = @_;
	$self -> {_log} = $log;	
}

sub raw2mzXML {
	my ($self, $mode) = @_;
	my $LOG = $self -> getLogFile();	
	my $rawFile = $self -> getRawFile();
	my $dirname = $self -> getRawFileDirname();
	my $rawFileBasename = $self -> getRawFileBasename();
	my $converter = $self -> getConverter();

	my $mzXML = "$dirname/$rawFileBasename.mzXML";
	if (-e $mzXML) {
		print "  mzXML file ($mzXML) already exists\n";
		print  $LOG "  mzXML file ($mzXML) already exists\n";		
	} else {
		print "  Converting the RAW file to mzXML format\n";
		print $LOG "  Converting the RAW file to mzXML format\n";		
		if ($mode eq "profile") {
			system (qq($converter "$rawFile" $mzXML >/dev/null 2>&1));
		} elsif ($mode eq "centroid") {
			system (qq($converter --mzXML -c "$rawFile" $mzXML >/dev/null 2>&1));		
		}
	}
	$self -> {_mzXMLFile} = $mzXML;
	return $mzXML;
}

sub getLogFile {
	my ($self) = @_;
	return $self -> {_log};	
}

sub getRawFile {
	my $self = shift;
	return $self -> {_rawFile};
}

sub getRawFileDirname {
	my $self = shift;
	my $rawFile = $self -> getRawFile();
	my @suffixlist = (".raw", ".RAW");
	my $rawFileDirname = dirname($rawFile, @suffixlist);
	return $rawFileDirname;
}

sub getRawFileBasename {
	my $self = shift;
	my $rawFile = $self -> getRawFile();
	my @suffixlist = (".raw", ".RAW");
	my $rawFileBasename = basename($rawFile, @suffixlist);
	return $rawFileBasename;
}

sub getConverter {
	my $self = shift;
	if (!defined($self -> {_converter})) {
		$self -> {_converter} = "wine /usr/local/bin/ReAdW.exe --mzXML";
	}
	return $self -> {_converter};
}

=head

sub set_converter
{
	my ($self,$converter)=@_;
	$self->{_converter}=$converter; 
}


sub get_mzXML_file
{
	my $self=shift;
	return $self->{_mzXML_file};
}

=cut

1;
