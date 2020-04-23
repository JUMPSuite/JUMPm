#!/usr/bin/perl

######### Job #################################################
#                                                             #
#       **************************************************    #  
#       **** Summary module for JUMPm   		      ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2015 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::Summary;

use strict;
use warnings;
use File::Basename;
use Spiders::XMLParser;
use Spiders::MathUtils;
use Chemistry::File::Formula;
use Storable;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.03;


@ISA	 = qw(Exporter);
@EXPORT = qw(set_log_file get_log_file generate_feature_table find_charge pair_summary missile_summary missile_summary_unlabel final_summary calculate_FDR unique_best_hit calculate_intensity_ratio octet_rule);
 

sub new{
    my ($class,%arg) = @_;
    my $self = {

    };
    bless $self, $class;
    return $self;
}
=head
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
=cut

sub set_log_file
{
	my ($self,$log)=@_;
	$self->{_log}=$log;	
}

sub get_log_file
{
	my ($self)=@_;
	return $self->{_log};	
}

sub generate_feature_table
{
	my ($self,$dta_path) = @_;
	my $xml = new Spiders::XMLParser();
	my $dta_path_mzXML = $dta_path;
	$dta_path_mzXML=~s/(.*)\.(\d+)/$1/;
	
	my $mzXML = $dta_path_mzXML . ".mzXML";	
	open (XML, "<$mzXML") || die "can not open the file";
    my $indexOffset = $xml->get_IndexOffset(*XML);
	
    my ($index_array, $last_scan) = $xml->get_IndexArray(*XML, $indexOffset);
	my %scan_RT;
    foreach my $index (@$index_array)
    {
		my ($scan) = $xml->get_scannum(*XML, $index);	
		my ($rt) = $xml->get_RT(*XML, $index);
		$scan_RT{$scan}=$rt;
	}
	
	open(FEATURE,"$dta_path.tmp.feature") || die "can not open the file:$!";
	open(FEATUREISO,">$dta_path.feature") || die "can not open the file:$!";
	print FEATUREISO "index","\t","m\/z","\t","z","\t","MS1 scan#","\t","RT","\t","Intensity","\t","S\/N","\n";
	my %hash;
	my %intensity_hash;
	<FEATURE>;
	while(<FEATURE>)
	{
		chomp $_;
		my @data=split(/\t/,$_);
		$hash{$data[2]}{$data[3]}{$data[1]}=$_;
		$intensity_hash{$data[1]}=$data[3];
	}
	close(FEATURE);	
	my %charge;
	my %isotope;
	foreach my $scan (keys %hash)
	{
		foreach my $intensity (reverse sort {$b<=>$a} keys %{$hash{$scan}})
		{
			foreach my $mz (keys %{$hash{$scan}{$intensity}})
			{	

				my $charge_hash = find_charge($mz,\%{$hash{$scan}});
				
				foreach my $charged_mz (keys %$charge_hash)
				{
					$charge{$scan}{$mz}=$charge_hash->{$charged_mz};
					if($intensity_hash{$mz}>$intensity_hash{$charged_mz})
					{
						$isotope{$scan}{$charged_mz}=$charge_hash->{$charged_mz};				
					}
					else
					{
						$isotope{$scan}{$mz}=$charge_hash->{$mz};
					}
				}			
			}
		}
	}
	my $index=0;
	foreach my $scan (sort {$a<=>$b} keys %hash)
	{
		foreach my $intensity (sort {$b<=>$a} keys %{$hash{$scan}})
		{
			foreach my $mz (sort {$a<=>$b} keys %{$hash{$scan}{$intensity}})
			{	
				next if(defined($isotope{$scan}{$mz}));
				my @data =split("\t",$hash{$scan}{$intensity}{$mz});
				$index++;
				if(defined($charge{$scan}{$mz}))
				{
					print FEATUREISO $index,"\t",sprintf("%.4f",$mz),"\t",$charge{$scan}{$mz},"\t",$scan,"\t",$scan_RT{$scan},"\t",sprintf("%.0f",$intensity),"\t",sprintf("%.1f",$data[4]),"\n";			
				}
				else
				{
					print FEATUREISO $index,"\t",sprintf("%.4f",$mz),"\t","0","\t",$scan,"\t",$scan_RT{$scan},"\t",sprintf("%.0f",$intensity),"\t",sprintf("%.1f",$data[4]),"\n";						
				}	
			}
		}
	}	
	close(FEATUREISO);
}

sub find_charge
{

	my ($select_mz,$hash)=@_;
	my $C=1.00335;
	#my $parameter = get_parameter();
	#my $tolerance = $parameter->{'decharge_ppm'};
	my $intrappm=10;
	my $maxcharge=6;
	my ($previous_peak_mz_low,$previous_peak_mz_high) =  ($select_mz, $select_mz + $C); 
	my %charge_hash;	
	foreach my $intensity (keys %$hash)
	{

		foreach my $mz ( keys %{$hash->{$intensity}})
		{

		# search the previous peak ( only one peak) 
			if($mz>$previous_peak_mz_low && $mz<$previous_peak_mz_high)
			{
				
				my $diff = 1/abs($mz-$select_mz);
				my $round_diff = sprintf("%.0f", $diff);
				next if ($round_diff==0);
				next if ($round_diff > $maxcharge);
				my $var = abs(abs($mz-$select_mz)-($C/$round_diff));
				next if ($var > ($intrappm/1000000)*$select_mz);
		                
				my $diffsum += ($var*1000000)/$mz;
				my $charge = $round_diff;
				$charge_hash{$mz}=$charge;
			}
		}
	}
	return \%charge_hash;
}

sub pair_summary
{
	my ($self,$dta_path) = @_;
	my $LOG = $self->get_log_file();
	my %unique_pair;
	my $pair_file = $dta_path . ".pair.all.comb";
	my $pair_file_uniq = $dta_path . ".pair.all.comb.uniq";
	
	open(PAIR,">$pair_file");
	print PAIR "Scan\tC12 mass\tN15 mass\tC13 mass\tC number\tN number\tPair score\n";
	open(UNIQPAIR,">$pair_file_uniq");
	print UNIQPAIR "Scan\tC12 mass\tN15 mass\tC13 mass\tC number\tN number\tPair score\n";
	
	my @pairs = glob("$dta_path/*.pair.all.comb");	
	foreach (@pairs)
	{
		open(FILE,$_);
		my $filename =  basename($_,  ".MS1.pair.all.comb");
		
		while(<FILE>)
		{
			chomp $_;
			my @line = split(/\t/,$_);
			if($#line>3)
			{
				my @C_num = split(/\:/,$line[1]);
				my @N_num = split(/\:/,$line[2]);
				print PAIR $filename,"\t",$line[0],"\t",$N_num[1],"\t",$C_num[1],"\t",$C_num[0],"\t",$N_num[0],"\t",$line[4],"\n";
				my $pair_line = $filename . "\t" . $line[0] . "\t" . $N_num[1] . "\t" . $C_num[1] . "\t" . $C_num[0] . "\t" . $N_num[1] . "\t" . $line[4] . "\n";
				$unique_pair{$line[0] . $N_num[1] . $C_num[1]}=$pair_line;	
			}
		}		
	}
	foreach (keys %unique_pair)
	{
		print UNIQPAIR $unique_pair{$_};
	}
	close(PAIR);
	close(UNIQPAIR);
	my $number_pair_CNC = scalar keys %unique_pair;
	print "  $number_pair_CNC three-peak pair(s) (C-N-C) were detected\n";
	print $LOG "  $number_pair_CNC three-peak pair(s) (C-N-C) were detected\n";	
}

sub missile_summary
{
	my ($self,$dta_path) = @_;
	my $LOG = $self->get_log_file();
	my $ms_hash_mol;

	my $missile_file = $dta_path . ".missile";
	open(MISSILE,">$missile_file");
	print MISSILE "Scan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tN15 Intensity\tC13 Intensity\tFormula\tAdduct\n";
	open(UNIQUEMISSILE,">$dta_path.formula");
	print UNIQUEMISSILE "Index\tScan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tN15 Intensity\tC13 Intensity\tFormula\tAdduct\tType\n";
	open(UNIQUEMISSILEPASS,">$missile_file.unique.pass");
	print UNIQUEMISSILEPASS "Scan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tN15 Intensity\tC13 Intensity\tFormula\tAdduct\tType\n";	

#	print "\n  Generating a summary table containing peak pairing information (missiles) information\n";
#	print $LOG "\n  Generating a summary table containing peak pairing information (missiles) information\n";	
	my @missile = glob("$dta_path/*.missile");
	my %uniq_missile;
	my %normal_missile;
	my %decoy_missile;
	my %uniq_missile_flag;
	my %FDR_distribution;
	my %FDR_database_distribution;
	
	my %normal_missile_flag;
	my %decoy_missile_flag;	
	my %target_formula;
	my %decoy_formula;
	my $missile_count = 0;
	my $normal_count = 0;
	my $decoy_count = 0;
	my $normal_count_CC_pair = 0;
	my $decoy_count_CC_pair = 0;
	my $missile_all_target_CNC=0;
	my $missile_all_target_CC=0;
	my $missile_all_decoy_CNC = 0;
	my $missile_all_decoy_CC = 0;
	my $formula_table_index=0;
	my $adduct = "";
	foreach (@missile)
	{
		open(FILE,$_);
		<FILE>;
		while(<FILE>)
		{
			chomp $_;
			my @data = split(/\t/,$_);
			for(my $i=0;$i<=$#data;$i++)
			{
				if($data[$i] eq "" or $data[$i] eq "NA")
				{
					$data[$i] = 0;
				}
			}			
			#next if(!defined($data[16]));
			
			#$data[16]=~s/(\D+)1(\D+)/$1$2/g;
			#$data[16]=~s/(\D+)1$/$1/;
			if(defined($data[17]))
			{
				$adduct = $data[17];
			}
			else
			{
				$adduct = "";
			}
			my $scan = $data[11];
			
			$scan =~s/.*\/(\d+)\.MS1/$1/;
			$scan =~s/\.iso//;
			$scan =~s/\.MS1//;	
			
			my $key = join(":",@data[1..3]);
			my $score = $data[4];

			
			print MISSILE $scan,"\t",join("\t",@data[1..4]),"\t",$data[12],"\t",$data[13],"\t",$data[14],"\t",$data[16],"\t",$adduct,"\n";
			
			if($data[15] eq "pass")
			{
				my $C12_mass = sprintf("%.3f",$data[1]);

				my $N15_mass = 0;

				if($data[2] ne "NA")
				{
					$N15_mass = sprintf("%.3f",$data[2]);

				}
				my $C13_mass = 0;
				if($data[3] ne "NA")
				{
					$C13_mass = sprintf("%.3f",$data[3]);				
				}
				next if($uniq_missile_flag{$C12_mass}{$N15_mass}{$C13_mass}{$data[16]});
				
				$missile_count++;


				my $target_decoy = $self->octet_rule($data[16]);
				$uniq_missile_flag{$C12_mass}{$N15_mass}{$C13_mass}{$data[16]} = $target_decoy;
				$FDR_distribution{$data[4]}{$target_decoy}++;
				
				print UNIQUEMISSILEPASS $scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\t","pass.",$target_decoy,"\n";				
				if($target_decoy eq "Target")
				{
					if($N15_mass!=0)
					{
						$missile_all_target_CNC++;
					}
					else
					{
						$missile_all_target_CC++;						
					}
				}
				elsif($target_decoy eq "Decoy")
				{
					if($N15_mass!=0)
					{
						$missile_all_decoy_CNC++;
					}
					else
					{
						$missile_all_decoy_CC++;						
					}
				}				
			}
			
			if($data[15] eq "Target")
			{
				my $C12_mass = sprintf("%.3f",$data[1]);
				my $N15_mass = 0;
				if($data[2] ne "NA")
				{
					$N15_mass = sprintf("%.3f",$data[2]);
				}
				my $C13_mass = 0;
				if($data[3] ne "NA")
				{
					$C13_mass = sprintf("%.3f",$data[3]);				
				}
				
				next if($normal_missile_flag{$C12_mass}{$N15_mass}{$C13_mass}{$data[16]});
				$normal_count++;
				$formula_table_index++;
				print UNIQUEMISSILE $formula_table_index,"\t",$scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\t",$data[15],"\n";
				print UNIQUEMISSILEPASS $scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\tDatabase.",$data[15],"\n";					
				$normal_missile_flag{$C12_mass}{$N15_mass}{$C13_mass}{$data[16]} = $data[15];
				$FDR_database_distribution{$data[4]}{$data[15]}++;				
				$target_formula{$data[16]}=1;
				if($N15_mass==0)
				{
					$normal_count_CC_pair++;
				}				
			}
			elsif($data[15] eq "Decoy")
			{
				my $C12_mass = sprintf("%.3f",$data[1]);
				my $N15_mass = 0;
				if($data[2] ne "NA")
				{
					$N15_mass = sprintf("%.3f",$data[2]);
				}
				my $C13_mass = 0;
				if($data[3] ne "NA")
				{
					$C13_mass = sprintf("%.3f",$data[3]);				
				}
				
				next if($decoy_missile_flag{$C12_mass}{$N15_mass}{$C13_mass}{$data[16]});
				$decoy_count++;
				$formula_table_index++;
				print UNIQUEMISSILE $formula_table_index,"\t",$scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\t",$data[15],"\n";
				print UNIQUEMISSILEPASS $scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\tDatabase.",$data[15],"\n";				
				$decoy_missile_flag{$C12_mass}{$N15_mass}{$C13_mass}{$data[16]} = $data[15];
				$FDR_database_distribution{$data[4]}{$data[15]}++;					
				$decoy_formula{$data[16]}=1;
				if($N15_mass==0)
				{
					$decoy_count_CC_pair++;
				}				
			}
			
		}
		close(FILE);
	}

	my $MISSILE_FDR = 0;
	if($normal_count>0)
	{
		$MISSILE_FDR = sprintf("%.2f",$decoy_count / $normal_count * 100);
	}
	my $target_formula_number = scalar keys %target_formula;
	my $decoy_formula_number = scalar keys %decoy_formula;
	print "\n  Based on the searching database:";
	print $LOG "\n  Based on the searching database:";	
	print "\n  FDR estimation using defined database (formula level): $normal_count targets and $decoy_count decoys were found, with an FDR of ${MISSILE_FDR}%\n";	
	print $LOG "\n  FDR estimation using defined database (formula level): $normal_count targets and $decoy_count decoys were found, with an FDR of ${MISSILE_FDR}%\n";
	
	my $normal_count_CNC_pair = $normal_count - $normal_count_CC_pair;
	my $decoy_count_CNC_pair = $decoy_count - $decoy_count_CC_pair;		
	print "    Among $normal_count target hits: $normal_count_CNC_pair with three peak pairs (C-N-C), $normal_count_CC_pair with two peak pairs (C-C)\n";
	print "    Among $decoy_count decoy hits: $decoy_count_CNC_pair with three peak pairs (C-N-C), $decoy_count_CC_pair with two peak pairs (C-C)\n";
	
	print $LOG "    Among $normal_count target hits: $normal_count_CNC_pair with three peak pairs (C-N-C), $normal_count_CC_pair with two peak pairs (C-C)\n";
	print $LOG "    Among $decoy_count decoy hits: $decoy_count_CNC_pair with three peak pairs (C-N-C), $decoy_count_CC_pair with two peak pairs (C-C)\n";
		
	my $MISSILE_FDR_CNC = 0.00;
	my $MISSILE_FDR_CC = 0.00;
	if($normal_count_CNC_pair>0)
	{
		$MISSILE_FDR_CNC = sprintf("%.2f",$decoy_count_CNC_pair / $normal_count_CNC_pair * 100);
	}
	if($normal_count_CC_pair>0)
	{
		$MISSILE_FDR_CC = sprintf("%.2f",$decoy_count_CC_pair / $normal_count_CC_pair * 100);		
	}
	print "    FDR for three-peak pair (C-N-C) is ${MISSILE_FDR_CNC}%, and for two-peak pair (C-C) is ${MISSILE_FDR_CC}%, respectively\n";	
	print $LOG "    FDR for three-peak pair (C-N-C) is ${MISSILE_FDR_CNC}%, and for two-peak pair (C-C) is ${MISSILE_FDR_CC}%, respectively\n";		
	$target_formula_number = 1 if ($target_formula_number==0);
	my $formula_FDR = sprintf("%.2f",($decoy_formula_number / $target_formula_number) * 100);	
	$target_formula_number = 0 if ($target_formula_number==1);	
	print "  $target_formula_number target formula(s) and $decoy_formula_number decoy formula(s) with an FDR of ${formula_FDR}%\n";
	print  $LOG "  $target_formula_number target formula(s) and $decoy_formula_number decoy formula(s) with an FDR of ${formula_FDR}%\n";
	
	my ($sum_target,$sum_decoy,$cutoff_score) = $self->calculate_FDR(\%FDR_database_distribution,0.01);	
	print "  $sum_target target formula(s) and $sum_decoy decoy formula(s) with an FDR of 1%\n";
	print  $LOG "  $sum_target target formula(s) and $sum_decoy decoy formula(s) with an FDR of 1% at Pscore of $cutoff_score\n";
			
	my $missile_all_targets = $missile_all_target_CNC + $missile_all_target_CC + $normal_count;
	my $missile_all_decoys = $missile_all_decoy_CNC + $missile_all_decoy_CC + $decoy_count;	
	
	my $missile_all_FDR = sprintf("%.2f",($missile_all_decoy_CNC + $missile_all_decoy_CC) / ($missile_all_target_CNC + $missile_all_target_CC)*100);
=head	
	print "\n  Based on all hits passed simulated mass formula database, the FDR is ${missile_all_FDR}%\n";
	print $LOG "\n  Based on all hits passed simulated mass formula database, the FDR is ${missile_all_FDR}%\n";
	
	my $total_missile_pass_database = $missile_count + $normal_count + $decoy_count;	
	print "  A total of $total_missile_pass_database candidate MISSILES, including $missile_all_targets targets and $missile_all_decoys decoys\n";
	print $LOG "  A total of $total_missile_pass_database candidate MISSILES, including $missile_all_targets targets and $missile_all_decoys decoys\n";
	my $MISSILE_ALL_FDR_CNC = 0.00;
	my $MISSILE_ALL_FDR_CC = 0.00;
	$missile_all_target_CNC=$missile_all_target_CNC+$normal_count_CNC_pair;
	$missile_all_target_CC=$missile_all_target_CC+$normal_count_CC_pair;
	$missile_all_decoy_CNC=$missile_all_decoy_CNC+$decoy_count_CNC_pair;
	$missile_all_decoy_CC=$missile_all_decoy_CC+$decoy_count_CC_pair;
	
	if($missile_all_target_CNC>0)
	{
		$MISSILE_ALL_FDR_CNC = sprintf("%.2f",($missile_all_decoy_CNC / $missile_all_target_CNC) * 100);
	}
	if($missile_all_target_CC>0)
	{
		$MISSILE_ALL_FDR_CC = sprintf("%.2f",($missile_all_decoy_CC / $missile_all_target_CC) * 100);		
	}	

	print "    Among target hits: $missile_all_target_CNC with three-peak pairs (C-N-C), $missile_all_target_CC with two-peak pairs (C-C)\n";
	print "    Among decoy hits: $missile_all_decoy_CNC with three-peak pairs (C-N-C), $missile_all_decoy_CC with two-peak pairs (C-C)\n";
	print "    FDR for three-peak pair (C-N-C) is ${MISSILE_ALL_FDR_CNC}%, and for two-peak pair (C-C) is ${MISSILE_ALL_FDR_CC}%, respectively\n";	

		
	print $LOG "    Among target hits: $missile_all_target_CNC with three-peak pairs (C-N-C), $missile_all_target_CC with two-peak pairs (C-C)\n";
	print $LOG "    Among decoy hits: $missile_all_decoy_CNC with three-peak pairs (C-N-C), $missile_all_decoy_CC with two-peak pairs (C-C)\n";
	print $LOG "    FDR for three-peak pair (C-N-C) is ${MISSILE_ALL_FDR_CNC}%, and for two-peak pair (C-C) is ${MISSILE_ALL_FDR_CC}%, respectively\n";	
	($sum_target,$sum_decoy,$cutoff_score) = $self->calculate_FDR(\%FDR_distribution,0.01);	
	print "  $sum_target target formula(s) and $sum_decoy decoy formula(s) with an FDR of 1%\n";
	print  $LOG "  $sum_target target formula(s) and $sum_decoy decoy formula(s) with an FDR of 1%\n";	
=cut

	
	print "\n  Generating a summary table containing structural information\n";
	print $LOG "\n  Generating a summary table containing structural information\n";	
	
	my @temp_smout = glob("$dta_path/*.smout");
	my @smout;
	my %formula;
	my %smiles_structure;
	foreach my $smout_file (@temp_smout)
	{
		open(TFILE,$smout_file);
		my $cnt=0;
		my $head = <TFILE>;
		while(<TFILE>)
		{
			$cnt++;
			my @lines=split(/\t/,$_);
			$formula{$lines[9]}=1;
			$smiles_structure{$lines[12]}=1;			
		}
		next if($cnt<1);
		push (@smout,$smout_file);
	}
## format of smout
# Index	C12 mass	N15 mass	C13 mass	MS1 Scan	C12 Intensity	Formula	Pscore	Structure (SMILES) Type	
	my $num_scan_structure = scalar(@smout);
	my $uniq_formula_num = scalar keys (%formula);
	my $uniq_structure_num = scalar keys (%smiles_structure);
	print "  $uniq_formula_num formulas are found to have $uniq_structure_num structures (SMILES) in the structure database\n";
	print $LOG "  $uniq_formula_num formulas are found to have $uniq_structure_num structures (SMILES) in the structure database\n";
		
	my $num_structure = 0;	
	foreach (@smout)
	{
		print "  Reading file: $_ \r";
		open(FILE,$_);
		<FILE>;
		while(<FILE>)
		{
			
			chomp $_;
			my @data = split(/\t/,$_);
			my $scan = $data[4];
			$scan =~s/.*\/(\d+)\.MS1/$1/;
			$scan =~s/\.iso//;		
			$scan =~s/\.MS1//;
			for(my $i=0;$i<=$#data;$i++)
			{
				if($data[$i] eq "" or $data[$i] eq "NA")
				{
					$data[$i] = 0;
				}
			}			
			
			my $N15_mass = 0;
			if($data[2] ne "NA")
			{
				$N15_mass = sprintf("%.6f",$data[2]);				
			}
			my $C13_mass = 0;
			if($data[3] ne "NA")
			{
				$C13_mass = sprintf("%.6f",$data[3]);				
			}
			
			my $mass = sprintf("%.6f", $data[1]);
			
			$ms_hash_mol->{$scan}->{$mass}->{$data[12]}->{'C12'}=$_;
			$ms_hash_mol->{$scan}->{$N15_mass}->{$data[12]}->{'N15'}=$_;						
			$ms_hash_mol->{$scan}->{$C13_mass}->{$data[12]}->{'C13'}=$_;
			$num_structure++;	
		}
		close(FILE);
	}	
	
	return $ms_hash_mol;
}			


sub missile_summary_unlabel
{
	my ($self,$dta_path) = @_;
	my $LOG = $self->get_log_file();
	my $ms_hash_mol;

	my $missile_file = $dta_path . ".missile";
	open(MISSILE,">$missile_file");
	
=head	
	print MISSILE "Scan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tFormula\tAdduct\tType\n";
	open(UNIQUEMISSILE,">$missile_file.formula");
	print UNIQUEMISSILE "Index\tScan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tFormula\tAdduct\tType\n";
	open(UNIQUEMISSILEPASS,">$missile_file.unique.pass");
	print UNIQUEMISSILEPASS "Scan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tFormula\tAdduct\tType\n";	
=cut
	
	print MISSILE "Scan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tN15 Intensity\tC13 Intensity\tFormula\tAdduct\tType\n";
	open(UNIQUEMISSILE,">$dta_path.formula");
	print UNIQUEMISSILE "Index\tScan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tN15 Intensity\tC13 Intensity\tFormula\tAdduct\tType\n";
	open(UNIQUEMISSILEPASS,">$missile_file.unique.pass");
	print UNIQUEMISSILEPASS "Scan\tC12 mass\tN15 mass\tC13 mass\tPair score\tC12 Intensity\tN15 Intensity\tC13 Intensity\tFormula\tAdduct\tType\n";		
	
	
	print "\n  Generating a summary table containing structural information\n";
	my $formula_table_index=0;
	my @temp_smout = glob("$dta_path/*.smout");
	my @smout;
	my %formula;
	my %smiles_structure;
	foreach my $smout_file (@temp_smout)
	{
		open(TFILE,$smout_file);
		my $cnt=0;
		my $head = <TFILE>;
		while(<TFILE>)
		{
			$cnt++;
			my @lines=split(/\t/,$_);
			$formula{$lines[9]}=1;
			$smiles_structure{$lines[12]}=1;
			
		}
		next if($cnt<1);
		push (@smout,$smout_file);
	}
## format of smout
# Index	C12 mass	N15 mass	C13 mass	MS1 Scan	C12 Intensity	Formula	Pscore	Structure (SMILES) Type	
	my $num_scan_structure = scalar(@smout);
	my $uniq_formula_num = scalar keys (%formula);
	my $uniq_structure_num = scalar keys (%smiles_structure);
	print "  A total of $num_scan_structure MS1 scans, $uniq_formula_num unique formulas, $uniq_structure_num unique structures were found\n";
	print $LOG "  A total of $num_scan_structure MS1 scans, $uniq_formula_num unique formulas, $uniq_structure_num unique structures were found\n";
	
	my $num_structure = 0;	
	foreach (@smout)
	{
		print "  Reading file: $_ \r";
		open(FILE,$_);
		<FILE>;

		while(<FILE>)
		{
			
			chomp $_;
			my @data = split(/\t/,$_);
			my $scan = $data[4];
			if($scan =~/.*\/(\d+)\.MS1/)
			{
				$scan =~s/.*\/(\d+)\.MS1/$1/;
			}
			$scan =~s/\.iso//;		
			$scan =~s/\.MS1//;

			for(my $i=0;$i<=$#data;$i++)
			{
				if($data[$i] eq "" or $data[$i] eq "NA")
				{
					$data[$i] = 0;
				}
			}			
					
			my $mass = sprintf("%.6f", $data[1]);
			
			$ms_hash_mol->{$scan}->{$mass}->{$data[12]}->{'C12'}=$_;
			$num_structure++;	
		}
		close(FILE);
	}

	print "\n\n  Generating a summary table containing formula information\n";
	print $LOG "\n\n  Generating a summary table containing formula information\n";	
	my @missile = glob("$dta_path/*.missile");
	my %uniq_missile;
	my %normal_missile;
	my %decoy_missile;
	my %uniq_missile_flag;
	my %normal_missile_flag;
	my %decoy_missile_flag;	
	my $missile_count = 0;
	my $normal_count = 0;
	my $decoy_count = 0;
	my $adduct="";
	foreach (@missile)
	{
		open(FILE,$_);
		<FILE>;
		while(<FILE>)
		{
			chomp $_;
			my @data = split(/\t/,$_);
			my $scan = $data[11];
			
			$scan =~s/.*\/(\d+)\.MS1/$1/;
			$scan =~s/\.iso//;
			$scan =~s/\.MS1//;	
			for(my $i=0;$i<=$#data;$i++)
			{
				if($data[$i] eq "" or $data[$i] eq "NA")
				{
					$data[$i] = 0;
				}
			}			
			if(defined($data[17]))
			{
				$adduct = $data[17];
			}
			else
			{
				$adduct = "";
			}
			$data[16]=~s/(\D+)1(\D+)/$1$2/g;
			$data[16]=~s/(\D+)1$/$1/;
			
			my $key = join(":",@data[1..3]);
			my $score = $data[4];
=head
			my $target_decoy = $self->octet_rule($data[16]);			
			#print MISSILE $scan,"\t",join("\t",@data[1..4]),"\t",$data[12],"\t",$data[14],"\t",$data[13],"\n";
			print UNIQUEMISSILEPASS $scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$adduct,"\t",$data[15],"\t","pass.",$target_decoy,"\n";				
			if($data[15] eq "pass")
			{
				my $C12_mass = sprintf("%.3f",$data[1]);
		
				next if($uniq_missile_flag{$C12_mass}{$data[16]});
				
				$missile_count++;
				print UNIQUEMISSILEPASS $scan,"\t",join("\t",@data[1..4]),"\t",$data[12],"\t",$data[14],"\t",$data[13],"\n";
				
				$uniq_missile_flag{$C12_mass}{$data[16]} = 1;
			}
=cut			
			if($data[15] eq "Target")
			{
				my $C12_mass = sprintf("%.3f",$data[1]);
				my $N15_mass = 0;
				if($data[2] ne "NA")
				{
					$N15_mass = sprintf("%.3f",$data[2]);				
				}
				my $C13_mass = 0;
				if($data[3] ne "NA")
				{
					$C13_mass = sprintf("%.3f",$data[3]);				
				}
				if($data[13] eq "NA")
				{
					$data[13] = 0;				
				}
				if($data[14] eq "NA")
				{
					$data[14] = 0;				
				}
				
				next if($normal_missile_flag{$C12_mass}{$data[16]});
				$normal_count++;
				$formula_table_index++;
				#print UNIQUEMISSILEPASS $scan,"\t",join("\t",@data[1..4]),"\t",$data[12],"\t",$data[14],"\t",$data[13],"\n";				
				#print UNIQUEMISSILE $formula_table_index,"\t",$scan,"\t",join("\t",@data[1..4]),"\t",$data[12],"\t",$data[14],"\t",$data[13],"\n";	

				print UNIQUEMISSILE $formula_table_index,"\t",$scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\t",$data[15],"\n";
				print UNIQUEMISSILEPASS $scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\tpass.",$data[15],"\n";				
				
				$normal_missile_flag{$C12_mass}{$data[16]} = 1;
			}
			elsif($data[15] eq "Decoy")
			{
				my $C12_mass = sprintf("%.3f",$data[1]);
				my $N15_mass = 0;
				if($data[2] ne "NA")
				{
					$N15_mass = sprintf("%.3f",$data[2]);				
				}
				my $C13_mass = 0;
				if($data[3] ne "NA")
				{
					$C13_mass = sprintf("%.3f",$data[3]);				
				}
				if($data[13] eq "NA")
				{
					$data[13] = 0;				
				}
				if($data[14] eq "NA")
				{
					$data[14] = 0;				
				}
				
				next if($decoy_missile_flag{$C12_mass}{$data[16]});
				$decoy_count++;
				$formula_table_index++;
				#print UNIQUEMISSILEPASS $scan,"\t",join("\t",@data[1..4]),"\t",$data[12],"\t",$data[14],"\t",$data[13],"\n";				
				#print UNIQUEMISSILE $formula_table_index,"\t",$scan,"\t",join("\t",@data[1..4]),"\t",$data[12],"\t",$data[14],"\t",$data[13],"\n";

				print UNIQUEMISSILE $formula_table_index,"\t",$scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\t",$data[15],"\n";
				print UNIQUEMISSILEPASS $scan,"\t",sprintf("%.4f\t%.4f\t%.4f\t%.2f",@data[1..4]),"\t",sprintf("%.0f\t%.0f\t%.0f",@data[12..14]),"\t",$data[16],"\t",$adduct,"\tpass.",$data[15],"\n";				
				$decoy_missile_flag{$C12_mass}{$data[16]} = 1;
			}
			
		}
		close(FILE);
	}

	my $MISSILE_FDR = 0;
	if($normal_count>0)
	{
		$MISSILE_FDR = sprintf("%.2f",$decoy_count / $normal_count * 100);
	}

	print "\n  FDR estimation using defined database (formula level): $normal_count targets and $decoy_count decoys were found, with an FDR of ${MISSILE_FDR}%\n";	
	
	print $LOG "\n  FDR estimation using defined database (formula level): $normal_count targets and $decoy_count decoys were found, with an FDR of ${MISSILE_FDR}%\n";	

	return $ms_hash_mol;
}			


sub final_summary
{
	my ($self,$ms_hash_mol,$MS1_MS2_matched,$dta_path) = @_;
	my $missile_structure_with_score = 0;
	my $missile_structure_without_score = 0;
	
	my $LOG = $self->get_log_file();
	my $MS2_structure = $dta_path . ".MS2" . ".structure";
######### output file ##################################
	open(OUTPUT,">$dta_path.spectrum_matches") || die "can not open the summary file"; 	
	open(OUTPUTBEST,">$MS2_structure") || die "can not open the summary file";
	print OUTPUT "Index\tMS1 scan\tC12 mass\tN15 mass\tC13 mass\tC12 intensity\tN15 intensity\tC13 intensity\tPscore\tFormula\tAdduct\tTarget-Decoy\tMS2 scan\tName\tSMILES\tInChIKey\tMscore\tPrecursor Type\n";
	print OUTPUTBEST "MS1 scan\tC12 mass\tN15 mass\tC13 mass\tC12 intensity\tN15 intensity\tC13 intensity\tPscore\tFormula\tAdduct\tTarget-Decoy\tMS2 scan\tName\tSMILES\tInChIKey\tMscore\tPrecursor Type\n";
######### programming starting information #############
	my $index=0;
	my $best_index=0;
	print "\n\n  Generating final summary results\n";
	print $LOG "\n\n  Generating final summary results\n";	
	foreach my $scan_missile (sort {$a<=>$b} keys %$ms_hash_mol)
	{

		foreach my $mz_missile (keys %{$ms_hash_mol->{$scan_missile}})
		{
			my $besthit = 1;	
			next if ($mz_missile == 0);
			my $flag = 0;
			my %score_hash;		

			foreach my $MS2_scan (keys %{$MS1_MS2_matched->{$scan_missile}->{$mz_missile}})
			{

				#my $MS2_scan = $MS1_MS2_matched->{$scan_missile}->{$mz_missile};
				my $score_file = $MS2_scan . ".MS2.score";
				print "  summarizing scoring file: $score_file \r";
				next if(!-e("$dta_path/$score_file"));
				$flag = 1;						
				open(SCORE,"$dta_path/$score_file") || die "  can not open the file: $score_file \n";

###### format of score file
#MS2 scan	SMILES	Matched # peaks	Total theoretical #peak	MScore	Precursor type	Rank score	Intensity ratio
				
				while(<SCORE>)
				{
					chomp $_;			
					my @lines = split(/\t/,$_);
					next if(!defined($lines[1]));					
					$lines[1]=~s/\\//g;
#					my $MS2_scan = basename($lines[0]);
					next if(!defined($lines[4]));	
					$score_hash{$lines[4]}{'MS2_scan'} = $MS2_scan ;
					push (@{$score_hash{$lines[4]}{'smiles'}}, $_);					
				}
				close(SCORE);
			}
# pvalue is score				
			foreach my $pvalue (sort {$b<=>$a} keys %score_hash)
			{
				foreach (@{$score_hash{$pvalue}{'smiles'}})
				{
###### format of smout file
#NORM	1	120.080260605046	121.077216674583	128.106971970358	/home/xwang4/JUMPm/comparison_noFClBr/labeled_10ppm/yeast_ms2_score/yeast_ms2_score.54/1002.MS1	141843694.428125	C8H9N1		C1C2CC(C1C=C2)C#N C8H9N					
					my $score_smile = $_;
					my @score_line=split(/\t/,$score_smile);
					my $ms2_scan = $score_line[0];			
					my $mscore= $score_line[4];
					my $prec_type = $score_line[5];
					my $rank_score= $score_line[6];
					
######### a bug fixed in version 1.6.6: smiles also contains "\" 

					$score_line[1]=~s/\\\(/\(/g;
					$score_line[1]=~s/\\\)/\)/g;

					my $smout = $ms_hash_mol->{$scan_missile}->{$mz_missile}->{$score_line[1]};				

					my ($smout_value) = values %{$smout};
					next if(!defined($smout_value));
					my @smout_line = split(/\t/,$smout_value);
										
					if (!defined($smout_line[4]))
					{
							next;
					}
					my $scan = basename($smout_line[4]);
					my $N15_intensity = 0;
					my $C13_intensity = 0;
					if($smout_line[2]!=0)
					{
						$N15_intensity =$smout_line[6];
					}
					if(defined($smout_line[7]))
					{
						$C13_intensity =$smout_line[7];
					}						
					$scan =~s/(\d+)\.MS1/$1/;
					$scan =~s/\.iso//;	
					$index++;
					print OUTPUT $index,"\t",$scan,"\t",sprintf("%.4f",$smout_line[1]),"\t",sprintf("%.4f",$smout_line[2]),"\t",sprintf("%.4f",$smout_line[3]),"\t",$smout_line[5],"\t",sprintf("%.0f",$N15_intensity),"\t",sprintf("%.0f",$smout_line[7]),"\t",sprintf("%.2f",$smout_line[8]),"\t",$smout_line[9],"\t",$smout_line[10],"\t",$smout_line[14],"\t",basename($score_hash{$pvalue}{'MS2_scan'}),"\t",$smout_line[11],"\t",$smout_line[12],"\t",$smout_line[13],"\t",sprintf("%.2f",$pvalue),"\t",$prec_type,"\n";	
					
					
					if($besthit==1)
					{
						$missile_structure_with_score++;
						$best_index++;
						print OUTPUTBEST $scan,"\t",sprintf("%.4f",$smout_line[1]),"\t",sprintf("%.4f",$smout_line[2]),"\t",sprintf("%.4f",$smout_line[3]),"\t",$smout_line[5],"\t",sprintf("%.0f",$N15_intensity),"\t",sprintf("%.0f",$smout_line[7]),"\t",sprintf("%.2f",$smout_line[8]),"\t",$smout_line[9],"\t",$smout_line[10],"\t",$smout_line[14],"\t",basename($score_hash{$pvalue}{'MS2_scan'}),"\t",$smout_line[11],"\t",$smout_line[12],"\t",$smout_line[13],"\t",sprintf("%.2f",$pvalue),"\t",$prec_type,"\n";

						$besthit=0;						
					}						
				}
			}				
		}
	}
	
	close(OUTPUT);
	close(OUTPUTBEST);	
	
	print "  $missile_structure_with_score unique MISSILES matched to MS2 scans\n";
	print $LOG "  $missile_structure_with_score unique MISSILES matched to MS2 scans\n";
	
	print "\n";
	$self->unique_best_hit($MS2_structure);
}

sub final_summary_dta
{
	my ($self,$ms_hash_mol,$MS1_MS2_matched,$dta_path) = @_;
	my $missile_structure_with_score = 0;
	my $missile_structure_without_score = 0;
	
	my $LOG = $self->get_log_file();
	my $MS2_structure = $dta_path . ".MS2" . ".structure";
######### output file ##################################
	open(OUTPUT,">$dta_path.spectrum_matches") || die "can not open the summary file"; 	
	open(OUTPUTBEST,">$MS2_structure") || die "can not open the summary file";
	print OUTPUT "Index\tMS1 scan\tC12 mass\tN15 mass\tC13 mass\tC12 intensity\tN15 intensity\tC13 intensity\tPscore\tFormula\tAdduct\tTarget-Decoy\tMS2 scan\tName\tSMILES\tInChIKey\tMscore\tPrecursor Type\n";
	print OUTPUTBEST "MS1 scan\tC12 mass\tN15 mass\tC13 mass\tC12 intensity\tN15 intensity\tC13 intensity\tPscore\tFormula\tAdduct\tTarget-Decoy\tMS2 scan\tName\tSMILES\tInChIKey\tMscore\tPrecursor Type\n";
######### programming starting information #############
	my $index=0;
	my $best_index=0;
	print "\n\n  Generating final summary results\n";
	print $LOG "\n\n  Generating final summary results\n";	
	foreach my $scan_missile (keys %$ms_hash_mol)
	{

		foreach my $mz_missile (keys %{$ms_hash_mol->{$scan_missile}})
		{
			my $besthit = 1;	
			next if ($mz_missile == 0);
			my $flag = 0;
			my %score_hash;		

			foreach my $MS2_scan (keys %{$MS1_MS2_matched->{$scan_missile}->{$mz_missile}})
			{

				#my $MS2_scan = $MS1_MS2_matched->{$scan_missile}->{$mz_missile};
				my $score_file = $MS2_scan . ".score";
				print "  summarizing scoring file: $score_file \r";
				next if(!-e("$dta_path/$score_file"));
				$flag = 1;						
				open(SCORE,"$dta_path/$score_file") || die "  can not open the file: $score_file \n";

###### format of score file
#MS2 scan	SMILES	Matched # peaks	Total theoretical #peak	MScore	Precursor type	Rank score	Intensity ratio
				
				while(<SCORE>)
				{
					chomp $_;			
					my @lines = split(/\t/,$_);
					next if(!defined($lines[1]));					
					$lines[1]=~s/\\//g;
#					my $MS2_scan = basename($lines[0]);
					next if(!defined($lines[4]));	
					$score_hash{$lines[4]}{'MS2_scan'} = $MS2_scan ;
					push (@{$score_hash{$lines[4]}{'smiles'}}, $_);					
				}
				close(SCORE);
			}
# pvalue is score				
			foreach my $pvalue (sort {$b<=>$a} keys %score_hash)
			{
				foreach (@{$score_hash{$pvalue}{'smiles'}})
				{
###### format of smout file
#NORM	1	120.080260605046	121.077216674583	128.106971970358	/home/xwang4/JUMPm/comparison_noFClBr/labeled_10ppm/yeast_ms2_score/yeast_ms2_score.54/1002.MS1	141843694.428125	C8H9N1		C1C2CC(C1C=C2)C#N C8H9N					
					my $score_smile = $_;
					my @score_line=split(/\t/,$score_smile);
					my $ms2_scan = $score_line[0];			
					my $mscore= $score_line[4];
					my $prec_type = $score_line[5];
					my $rank_score= $score_line[6];
					
######### a bug fixed in version 1.6.6: smiles also contains "\" 

					$score_line[1]=~s/\\\(/\(/g;
					$score_line[1]=~s/\\\)/\)/g;

					my $smout = $ms_hash_mol->{$scan_missile}->{$mz_missile}->{$score_line[1]};				

					my ($smout_value) = values %{$smout};
					next if(!defined($smout_value));
					my @smout_line = split(/\t/,$smout_value);
										
					if (!defined($smout_line[4]))
					{
							next;
					}
					my $scan = basename($smout_line[4]);
					my $N15_intensity = 0;
					my $C13_intensity = 0;
					if($smout_line[2]!=0)
					{
						$N15_intensity =$smout_line[6];
					}
					if(defined($smout_line[7]))
					{
						$C13_intensity =$smout_line[7];
					}						

					$scan =~s/\.iso//;	
					$index++;
					print OUTPUT $index,"\t",$scan,"\t",sprintf("%.4f",$smout_line[1]),"\t",sprintf("%.4f",$smout_line[2]),"\t",sprintf("%.4f",$smout_line[3]),"\t",$smout_line[5],"\t",sprintf("%.0f",$N15_intensity),"\t",sprintf("%.0f",$smout_line[7]),"\t",sprintf("%.2f",$smout_line[8]),"\t",$smout_line[9],"\t",$smout_line[10],"\t",$smout_line[14],"\t",basename($score_hash{$pvalue}{'MS2_scan'}),"\t",$smout_line[11],"\t",$smout_line[12],"\t",$smout_line[13],"\t",sprintf("%.2f",$pvalue),"\t",$prec_type,"\n";	
					
					if($besthit==1)
					{
						$missile_structure_with_score++;
						$best_index++;
						print OUTPUTBEST $scan,"\t",sprintf("%.4f",$smout_line[1]),"\t",sprintf("%.4f",$smout_line[2]),"\t",sprintf("%.4f",$smout_line[3]),"\t",$smout_line[5],"\t",sprintf("%.0f",$N15_intensity),"\t",sprintf("%.0f",$smout_line[7]),"\t",sprintf("%.2f",$smout_line[8]),"\t",$smout_line[9],"\t",$smout_line[10],"\t",$smout_line[14],"\t",basename($score_hash{$pvalue}{'MS2_scan'}),"\t",$smout_line[11],"\t",$smout_line[12],"\t",$smout_line[13],"\t",sprintf("%.2f",$pvalue),"\t",$prec_type,"\n";

						$besthit=0;						
					}						
				}
			}				
		}
	}
	
	close(OUTPUT);
	close(OUTPUTBEST);	
	
	print "  $missile_structure_with_score unique MISSILES matched to MS2 scans\n";
	print $LOG "  $missile_structure_with_score unique MISSILES matched to MS2 scans\n";
	
	print "\n";
	$self->unique_best_hit($MS2_structure);
}

sub calculate_FDR
{
	my ($self,$FDR_distribution,$FDR_cutoff)=@_;
	
	my $cutoff_score=0;
	my $sum_target = 0;
	my $sum_decoy = 0;
	my $FDR = 0;
	foreach (reverse sort {$a<=>$b} keys %$FDR_distribution)
	{	
		#print $_,"\n";
		if($sum_target==0)
		{
			$FDR = 0;
		}
		else
		{
			$FDR = ($sum_decoy)/($sum_target);
		}
		if($FDR>$FDR_cutoff)
		{		
			last;			
		}
		if(defined($FDR_distribution->{$_}->{"Target"}))
		{
			$sum_target = $FDR_distribution->{$_}->{"Target"}+$sum_target;
		}
		if(defined($FDR_distribution->{$_}->{"Decoy"}))
		{
			$sum_decoy = $FDR_distribution->{$_}->{"Decoy"}+$sum_decoy;
		}
		$cutoff_score = $_;	

	}
	return ($sum_target,$sum_decoy,$cutoff_score);
}

sub unique_best_hit
{
	my ($self,$best_hit_file) = @_;
	my $LOG = $self->get_log_file();	
	my %best_hit_unique;
	my %score;	
	my $missile_structure_unique = 0;
	open(BESTHIT,$best_hit_file) || die "can not open the file: $!";
	
	$best_hit_file=~s/MS2\.structure//g;
	my $best_hit_unique_file = $best_hit_file . "structure"; 
	open(BESTHITUNIQ,">$best_hit_unique_file");
	my $head = <BESTHIT>;
	print BESTHITUNIQ "Index","\t",$head;
	while(<BESTHIT>)
	{
		chomp $_;
		my @line = split(/\t/,$_);
		next if($line[15] eq "NA");
		#change the missile into MS2 scan 
		#my $key = $line[1] . "_" . $line[2] . "_" . $line[3];
		my $key = $line[11];
		
		$best_hit_unique{$key}{$line[15]}= $_;						
	}
	close(BESTHIT);
	my $index=0;
	foreach my $key (keys %best_hit_unique)
	{
		foreach my $mscore (sort {$b<=>$a} keys %{$best_hit_unique{$key}})
		{
		
			$index++;
			print BESTHITUNIQ $index,"\t",$best_hit_unique{$key}{$mscore};
			$missile_structure_unique++;

			if($mscore<1)
			{
				$score{1}++;
			}
			elsif($mscore<5)
			{
				$score{5}++;								
			}
			elsif($mscore<10)
			{
				$score{10}++;								
			}
			else
			{
				$score{20}++;									
			}
			
			last;
		}
		print BESTHITUNIQ "\n";
	}
	close(BESTHITUNIQ);
	print "  $missile_structure_unique unique MISSILES with score\n";
	print $LOG "  $missile_structure_unique unique MISSILES with score\n";
	
	print "  MScore distribution\n";

	print  $LOG "  MScore distribution\n";

	
	foreach (sort {$a<=>$b} keys %score)
	{
		if($_==20)
		{
			print "  >=","10","\t",$score{$_},"\n";	
			print $LOG "  >=","10","\t",$score{$_},"\n";				
		}
		else
		{
			print "  <",$_,"\t",$score{$_},"\n";	
			print $LOG "  <",$_,"\t",$score{$_},"\n";
		}			
	}
	
}


sub calculate_intensity_ratio
{
	my ($self,$dta_path) = @_;
	my $LOG = $self->get_log_file();	
	my @missile = glob("$dta_path/*.pairs");
	my %ratio;
	my %defect;
	my @NC_int;
	my @CC_int;
	my @NC_defect;
	my @CC_defect;
	my $lib = $self->get_library_path();
	my $math = new Spiders::MathUtils();
	$math->set_library_path($lib);
	
	foreach (@missile)
	{
		open(FILE,$_);
		<FILE>;
		while(<FILE>)
		{
			chomp $_;
			my @data = split(/\t/,$_);

			
			if($data[2] != 0)
			{
				my $nc_defect = ($data[2]-$data[1])/int($data[2]-$data[1]+0.200000001);
				push(@NC_defect,$nc_defect);				
				push(@NC_int, log($data[5])/log(2));				
			}

			if($data[3] != 0)
			{			
				my $cc_defect = ($data[3]-$data[1])/int($data[3]-$data[1]+0.200000001);
				push(@CC_defect,$cc_defect);					
				push(@CC_int, log($data[7])/log(2));				
			}		
		}
		close(FILE);
	}
	if($#NC_defect<3)
	{
		
		print "\n\n  There is no enough isotopic peaks detected in this data set!\n";
		print "  Please make sure that your input is a lebeled dataset!\n\n";
		print $LOG "\n\n  There is no enough isotopic peaks detected in this data set!\n";
		print $LOG "  Please make sure that your input is a lebeled dataset!\n\n";		
		exit;
	}
	if($#NC_int<1)
	{
		$ratio{'NC'} = 0;
		$ratio{'CC'} = 0;
		$ratio{'NC_std'} = 0;
		$ratio{'CC_std'} = 0;
		$defect{'NC'}[0]=0;
		$defect{'CC'}[0]=0;
	
	}
	else
	{
		$ratio{'NC'} = $math->average(\@NC_int);
		$ratio{'CC'} = $math->average(\@CC_int);
		$ratio{'NC_std'} = $math->stdev(\@NC_int);
		$ratio{'CC_std'} = $math->stdev(\@CC_int);
		$defect{'NC'}[0] = $math->average(\@NC_defect);
		$defect{'CC'}[0] = $math->average(\@CC_defect);
		$defect{'NC'}[1] = $math->stdev(\@NC_defect);
		$defect{'CC'}[1] = $math->stdev(\@CC_defect);		
#		$defect{'NC'}=$math->est_cauchy_dist(\@NC_defect);
#		$defect{'CC'}=$math->est_cauchy_dist(\@CC_defect);
	}
	return (\%ratio,\%defect);
}

sub octet_rule
{
	my ($self,$f) = @_;
	my $target_decoy="";

	my %formula = Chemistry::File::Formula->parse_formula($f);
	if(!defined($formula{H}))
	{
		$formula{H}=0;
	}
	if(!defined($formula{Na}))
	{
		$formula{Na}=0;
	}
	if(!defined($formula{K}))
	{
		$formula{K}=0;
	}
	if(!defined($formula{N}))
	{
		$formula{N}=0;
	}
	if(!defined($formula{P}))
	{
		$formula{P}=0;
	}
	if(!defined($formula{F}))
	{
		$formula{F}=0;
	}
	if(!defined($formula{Cl}))
	{
		$formula{Cl}=0;
	}
	if(!defined($formula{Br}))
	{
		$formula{Br}=0;
	}
	if(!defined($formula{I}))
	{
		$formula{I}=0;
	}	
	my $hydrogen = $formula{H} + $formula{Na} + $formula{K} + $formula{N} + $formula{P} + $formula{F} + $formula{Cl} + $formula{Br} + $formula{I};
	if(($hydrogen % 2)==0)
	{
		$target_decoy = "Target";
	}
	else
	{
		$target_decoy = "Decoy";
	}
	return $target_decoy;
}

1;
