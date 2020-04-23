#!/usr/local/bin/perl 

use POSIX;
use File::Basename;
use Chemistry::File::Formula;
use Chemistry::MolecularMass;

my @input_dir = glob("/home/xwang4/JUMPm_database/MASS_FORMULA_DB/Formula{$ARGV[0]}*.txt");
#opendir(MASS_FORMULA_DB,$input_dir);
#while(my $file = readdir(MASS_FORMULA_DB))
foreach my $file (@input_dir)
{
	open(FILE,"$file");
	my $outfile=basename($file);
	
	open(OUT,">>new2/$outfile");
#	my $mass = $outfile;
#	$mass =~ s/Formula//;
	while(<FILE>)
	{
		chomp $_;
		my @data=split(/\t/,$_);
		if($data[1]=~/\_decoy/)
		{
		
			my @data1=split(/\:/,$data[0]);
			my $decoy_f = gen_decoy($data1[0]);
			my $mass = $data1[1]+2* 1.00782;
			my $mass_f = ceil($mass*100);
			open(OUT1,">>new2/Formula${mass_f}.txt");			
			print OUT1 $decoy_f,":",$mass,"\t",$data[1],"\n";
			close(OUT1);
		}
		else
		{
			print OUT $_,"\n";
		}
	}
	close(OUT);
}

sub gen_decoy
{
    my ($f)=shift;
	$f=~s/[+|-].*//;
	my %formula = Chemistry::File::Formula->parse_formula($f);
	my $hydrogen = $formula{H};
	my $new_hydrogen = $hydrogen + 2;
	$f =~ s/H$hydrogen/H$new_hydrogen/;
	return $f;
}

sub cal_mass
{
    my ($formula)=shift;

    my $mm = new Chemistry::MolecularMass;
    my $mass = sprintf("%.6f",$mm->calc_mass($formula));
	return $mass;
}



