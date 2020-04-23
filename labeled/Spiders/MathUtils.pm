#!/usr/bin/perl

## Release date: 01/31/2015
## Release version: version 11.1.1
## Module name: Spiders::MathUtils

######### Deisotope ##########################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################
# this value is used to remove the precursor ion

package Spiders::MathUtils;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);
use Statistics::R;

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT = qw(set_library_path get_library_path min  max  log10 combinations average  stdev est_cauchy_dist);

sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
}


sub set_library_path
{
	my ($self,$lib)=@_;
	$self->{_lib}=$lib;	
}

sub get_library_path
{
	my ($self)=@_;
	return $self->{_lib};	
}

sub min { 
	my ($self,$hash) = @_;

	my $min=1000000;
	foreach my $key (keys %$hash ) {
		$min = $hash->{$key} if $hash->{$key} < $min; 
	}
	return $min;
}

sub max { 
	my ($self,$hash) = @_;
	my $max=0;
	foreach my $key (keys %$hash ) { 
		$max = $hash->{$key} if $hash->{$key} > $max; 
	}
	return $max;
}

sub log10 
{
	my $n = shift;
	return log($n)/log(10);
}


sub combinations
{
    my ($self, $hash, $comb_num) = @_;
	
    my %new_hash;

    for (my $i = 0 ; $i < $comb_num ; $i++)
    {
        if ($i == 0)
        {
            foreach my $key (keys(%$hash))
            {
                $new_hash{$key} = $$hash{$key};
            }
        }
        else
        {
            foreach my $hash_key (keys(%new_hash))
            {
                foreach my $key (keys(%$hash))
                {
                    my $key_join = "$hash_key"."$key";
					my @alpha=sort(split(//,$key_join));
					my $key_comb=join('',@alpha);
                    $new_hash{$key_comb} = $$hash{$key} + $new_hash{$hash_key};
                }
            }
        }
    }
    return (%new_hash);
}

sub average {
	my ($self,$arrayref) = @_; # save the array passed to this function
	my $sum; # create a variable to hold the sum of the array's values
	foreach (@$arrayref) { $sum += $_; } # add each element of the array 
	# to the sum
	
	return $sum/@$arrayref; # divide sum by the number of elements in the
	# array to find the mean
}

sub stdev{
    my($self,$data) = @_;
    if(@$data == 1){
            return 0;
    }
    my $average = $self->average($data);
    my $sqtotal = 0;
    foreach(@$data) {
            $sqtotal += ($average-$_) ** 2;
    }
    my $std = ($sqtotal / (@$data-1)) ** 0.5;
    return $std;
}

sub est_cauchy_dist
{
	my ($self,$data) = @_;
	my $file = time();
	my $library = $self->get_library_path();
    my $R = Statistics::R->new(r_bin=>"$library/bin/R");

	open(TEMP,">/tmp/$file");
	print TEMP join("\n",@$data);
	print TEMP "\n";
	close(TEMP);
	$R->set('file',"/tmp/$file");
	$R->run(q`x<-read.table(file)`);

	$R->run(q`cauchy.fit <- function(theta,x){-sum(dcauchy(x,location=theta[1],scale=theta[2],log=TRUE),na.rm=T)}`,
	q`good <- !is.na(x)`,
	q`x = x[good]`,
	q`theta.start <- c(median(x),0.1)`,
	q`res <- nlminb(theta.start,cauchy.fit, x = x,lower=c(-10,1e-20),upper=c(10,10))`);
	my $result=$R->get('res$par');
	system(qq(rm -rf /tmp/$file));
	return $result;
}

1;

