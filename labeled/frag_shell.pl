#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin";
use Getopt::Long;
use Spiders::MS2_scoring;
use Spiders::Dta;
use Spiders::Params;
use Chemistry::File::Formula;
use Chemistry::File::SMILES;
use Chemistry::MolecularMass;


	
my ($help,$dtafile,$smile,$mass,$depth,$parameter,$prec_type);

GetOptions('-help|h'=>\$help,
		'-dtafile=s'=>\$dtafile,
		'-smile=s'=>\$smile,
		'-mass=s'=>\$mass,
		'-depth=s'=>\$depth,
		'-param=s'=>\$parameter,
		'-ptype=s'=>\$prec_type,
		'-target_decoy'=>\$target_decoy,
		);

my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();	

my $neutralLossFile = $Bin . "/neutralLoss.csv";	
my $bondEnergyFile =  $Bin . "/bondenergies.txt";
my $H = 1.007276466812;
my $dta=new Spiders::Dta;  
$dta->set_dta_file($dtafile) ;
$dta->process_dtafile($dtafile);
my $exp_mz = $dta->get_mz_array();

#my %mz_hash = %{$dta->{'_mz_hash'}};
my %mz_int_hash = %{$dta->{'_mz_int_hash'}};

#### MS2 peaks selection for comparison 
if(!defined($params->{'percentage_ms2_peaks'}))
{
	$params->{'percentage_ms2_peaks'} = 50;
}
my $percentage_MS2_peaks = $params->{'percentage_ms2_peaks'};
my $selected_peaks_number = int(($#$exp_mz+1) * $percentage_MS2_peaks/100);

my $number=0;
my @updated_mz;
my @updated_int;
my $sum_total_int = 0;

if($params->{'fragment_method'} == 0)
{
	unless(-e("$dtafile.score"))
	{
		open(SCORE,">$dtafile.score") || die "can not open the score file\n";
		print "MS2 file\tSMILES\tMatched peak number\tTheoretical peak number\tMscore\tPrecusor type\tRank_Mscore\tIntensity ratio\n";
		close(SCORE);
	}
	open(SCORE,">>$dtafile.score") || die "can not open the score file\n";
	flock SCORE,LOCK_EX;
	#print SCORE $dtafile,"\t",$smile,"\t", $matched,"\t",scalar(@frag_mz_values),"\t",$mscore,"\t",$prec_type,"\t",$wc_mscore,"\t",$int_ratio,"\n";
	print SCORE $dtafile,"\t",$smile,"\t", 0,"\t",0,"\t",0,"\t",$prec_type,"\t",0,"\t",0,"\n";

	close(SCORE);
	exit;
}


foreach my $mz (sort {$mz_int_hash{$b}<=>$mz_int_hash{$a}} keys %mz_int_hash)
{
	my $mz_int = int($mz);
	$mz_hash{$mz_int}{$mz} = $mz_int_hash{$mz};	
	$number++;
	push (@updated_mz,$mz);
	push (@updated_int,$mz_int_hash{$mz});
	$sum_total_int += $mz_int_hash{$mz}; 	
	last if($number>$selected_peaks_number);
}

my $scoring = new Spiders::MS2_scoring;
$scoring->set_exp_mz($exp_mz);
$scoring->set_parameter($params);
my @frag_mz_values=();
if($params->{'fragment_method'} == 1)
{
	## forth parameter is whether to break aromatic ring; 
	## fifth parameter is molecularFormulaRedundancyCheck; 
	##      * @param breakAromaticRings break aromatic rings?
	##     * @param molecularFormulaRedundancyCheck the experimental isomorphism check
	##  Before proceeding to the next fragment, a redundancy check is performed to eliminate duplicate fragments. 
	##  Redundancy occurs if a fragment A is part of both parent fragments AB and ABC, or the fragment A appears in different places of the molecule, as in ABA. 
	##  In both cases the redundant structures would cause longer runtimes and higher memory consumption without gaining any information. 
	##  In addition to full (and time consuming) graph isomorphism checks we describe simpler heuristics later in this paper.

	my @frag_return = SmileFragmentation($smile, $mass, $depth, "no", "yes", $neutralLossFile, $bondEnergyFile);
	exit if($frag_return[0]=~/Error/);

	my @tmp_frag_mz_values=();	
	for(my $i=2;$i<$#frag_return;$i++) 
	{
		next if ($frag_return[$i]=~/atom/);
		my @frag_data = split(/ /,$frag_return[$i]);
		print $frag_return[$i],"\t",$frag_data[1],"\t",$frag_data[-2],"\n";
		next if($frag_data[-2]=~/[CHNOPS]/);
		
		if($prec_type eq "N15")
		{
			my @data = split(/(C|N|H|O|P|S|F|Cl|Br)/,$frag_data[-3]);
			{
				my $pos_N = 0;	
				($pos_N) = grep {$data[$_] eq "N"} 0 .. $#data;
				my $num_N = 0; 			
				if($pos_N)
				{
					if($data[$pos_N+1]=="")
					{
						$num_N = 1;
					}
					else
					{
							$num_N = $data[$pos_N+1];
					}
				}
				$frag_data[-2] += 0.99703 * $num_N;
			}		
			push(@tmp_frag_mz_values,$frag_data[-2]);
		}
		elsif($prec_type eq "C13")
		{
			my @data = split(/(C|N|H|O|P|S|F|Cl|Br)/,$frag_data[-3]);
			{
				my $pos_C = 0;		
				($pos_C) = grep {$data[$_] eq "C"} 0 .. $#data;
				my $num_C = 0; 			
				if($pos_C)
				{
					if($data[$pos_C+1]=="")
					{
						$num_C = 1;
					}
					else
					{
						$num_C = $data[$pos_C+1];
					}
				}
				$frag_data[-2] += 1.00335 * $num_C;
			}		
			push(@tmp_frag_mz_values,$frag_data[-2]);
		}
		else
		{
			push(@tmp_frag_mz_values,$frag_data[-2]);	
		}	
	}

	if($params->{'mode'} == 1)
	{
		for(my $i=0;$i<=$#tmp_frag_mz_values;$i++)
		{
			#my $updated_frag_mz_values = $tmp_frag_mz_values[$i] + $H;
			my $updated_frag_mz_values = $tmp_frag_mz_values[$i] ;
			if($target_decoy eq "Decoy" and ($i % 2)==0)
			{
				$updated_frag_mz_values = $updated_frag_mz_values + $H;
			}
			push(@frag_mz_values,$updated_frag_mz_values)
		}
	}
	elsif($params->{'mode'} == -1)
	{
		for(my $i=0;$i<=$#tmp_frag_mz_values;$i++)
		{
			my $updated_frag_mz_values = $tmp_frag_mz_values[$i] - $H;
			if($target_decoy eq "Decoy" and ($i % 2)==0)
			{
				$updated_frag_mz_values = $updated_frag_mz_values - $H;
			}			
			push(@frag_mz_values,$updated_frag_mz_values)
		}
	}

}
elsif($params->{'fragment_method'} == 2)
{

	my $mm = new Chemistry::MolecularMass;
	my $command = "";
	$smile=~s/\\//g;
	my $outfile = time();
	if($params->{'mode'} == 1)
	{
		$command = "/usr/local/bin/wine $Bin/fraggraph-gen.exe \"$smile\" 1 + fragonly /tmp/$outfile 2>/dev/null";
	}
	else
	{
		$command = "/usr/local/bin/wine $Bin/fraggraph-gen.exe \"$smile\" 1 - fragonly /tmp/$outfile 2>/dev/null";	
	}
	system(qq($command));

	open(OUT,"/tmp/$outfile");
	while(<OUT>)
	{
        chomp $_;
        my @data =split(/\s+/, $_);
        my $mol = Chemistry::Mol->parse($data[2],format=>smiles);
        my $formula = $mol->print(format => formula);
		$mol = Chemistry::Mol->parse($formula,format=>"formula");		
		my $mass = $mol->mass;
        #my $mass = $mm->calc_mass($formula);
        next if($mass<50);

		if($prec_type eq "N15")
		{
			my @data = split(/(C|N|H|O|P|S|F|Cl|Br)/,$formula);
			{
				my $pos_N = 0;	
				($pos_N) = grep {$data[$_] eq "N"} 0 .. $#data;
				my $num_N = 0; 			
				if($pos_N)
				{
					if($data[$pos_N+1]=="")
					{
						$num_N = 1;
					}
					else
					{
						$num_N = $data[$pos_N+1];
					}
				}
				$mass += 0.99703 * $num_N;
			}		
			push(@frag_mz_values,$mass);
		}
		elsif($prec_type eq "C13")
		{
			my @data = split(/(C|N|H|O|P|S|F|Cl|Br)/,$formula);
			{
				my $pos_C = 0;		
				($pos_C) = grep {$data[$_] eq "C"} 0 .. $#data;
				my $num_C = 0; 			
				if($pos_C)
				{
					if($data[$pos_C+1]=="")
					{
						$num_C = 1;
					}
					else
					{
						$num_C = $data[$pos_C+1];
					}
				}
				$mass += 1.00335 * $num_C;
			}		
			push(@frag_mz_values,$mass);
		}
		else
		{
			push(@frag_mz_values,$mass);
		}
	}
	close(OUT);	
	system(qq(rm -rf /tmp/$outfile));
}
else
{
	print "please specify the right fragmentation method\n";
	exit;
}

my ($matched,$matched_hash_ref) = $scoring->compare_theoritical_experiment(\%mz_hash,\@frag_mz_values);
my $int_ratio = 0;
my $wc_p = 0;
my $mscore = 0;
my $wc_mscore = 0;




if($matched>0)
{
	my @sorted_updated_mz = sort {$a<=>$b} (@updated_mz);
	$scoring->set_exp_mz(\@sorted_updated_mz);
	my $p_value = $scoring->get_pvalue($matched,scalar(@frag_mz_values));

	## Get intensity for matched peaks
	my @matched_int_array;
	my $sum_matched_int = 0;
	foreach my $mz (keys %$matched_hash_ref)
	{
		push(@matched_int_array,$mz_int_hash{$mz});
		$sum_matched_int += $mz_int_hash{$mz};
	}
	

	
	$int_ratio = sprintf("%.2f",$sum_matched_int/$sum_total_int * 100);
	$wc_p = $scoring->wilcox_test(\@matched_int_array,\@updated_int);


	if($p_value != 0)
	{
		$mscore = -log($p_value) /log(10);
	}


	if($wc_p != 0)
	{
		$wc_mscore = -log($wc_p) /log(10);
	}
}

	
unless(-e("$dtafile.score"))
{
	open(SCORE,">$dtafile.score") || die "can not open the score file\n";
	print "MS2 file\tSMILES\tMatched peak number\tTheoretical peak number\tMscore\tPrecusor type\tRank_Mscore\tIntensity ratio\n";
	close(SCORE);
}
open(SCORE,">>$dtafile.score") || die "can not open the score file\n";
flock SCORE,LOCK_EX;
#print SCORE $dtafile,"\t",$smile,"\t", $matched,"\t",scalar(@frag_mz_values),"\t",$mscore,"\t",$prec_type,"\t",$wc_mscore,"\t",$int_ratio,"\n";
print SCORE $dtafile,"\t",$smile,"\t", $matched,"\t",scalar(@frag_mz_values),"\t",$mscore,"\t",$prec_type,"\t",join(";",keys (%$matched_hash_ref)),"\t",$int_ratio,"\n";

close(SCORE);


sub SmileFragmentation
{
	my ($smile, $mass, $depth, $breakAromaticRingFlag, $molecularFormulaRedundancyCheck, $neutralLossFile, $bondEnergyFile) = @_;
	if (!$smile) 
	{
		print "formula is not defined\n";
		return @array=();
	}
	if (!$mass) 
	{
		$mass = "0.0";
	}
	if (!$depth) 
	{
		$depth = "2";
	}
	if (!$molecularFormulaRedundancyCheck) 
	{
		$molecularFormulaRedundancyCheck = "no";
	}
	if (!$breakAromaticRingFlag) 
	{
			$breakAromaticRingFlag = "yes";
	}
	
	my $jarpath = getJarPath();
	#my $script = "java -jar " . $jarpath . " -FragmentSMILE " . $smile . " " . $mass . " " . $depth . " " . $breakAromaticRingFlag . " " . $molecularFormulaRedundancyCheck . " " . $neutralLossFile . " " . $bondEnergyFile;
	my $script = "java -jar " . $jarpath . " -AIMFragmenterWrapper " . "\"java -jar " . $jarpath . "\" .tmp " . $smile . " " . $mass . " " . $depth . " " . $breakAromaticRingFlag . " " . $molecularFormulaRedundancyCheck . " " . $neutralLossFile . " " . $bondEnergyFile;
	

	print $script;
    my @value = `$script`;
    return @value;
}

sub getJarPath
{
	my $JUMP_PATH = $Bin . "/JumpJarPackage";
	my $jarpath = $JUMP_PATH . "/JumpPackage.jar"; 
    return $jarpath;
}
