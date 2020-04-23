#!/usr/bin/perl  -I /home/jcho/dev/JUMPm/JUMPm_v1.8.0
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


#GetOptions('-help|h'=>\$help,
#		'-param=s'=>\$parameter,
#		'-peaks=s'=>\$peaks,
#		'-NC=f'=>\$NC,
#		'-CC=f'=>\$CC,
#		'-NC_std=f'=>\$NC_std,
#		'-CC_std=f'=>\$CC_std,
#		'-NC_defect_loc=f'=>\$NC_defect_loc,
#		'-CC_defect_loc=f'=>\$CC_defect_loc,
#		'-NC_defect_scale=f'=>\$NC_defect_scale,
#		'-CC_defect_scale=f'=>\$CC_defect_scale,
#		);
		
#my $p = Spiders::Params->new('-path'=>$parameter);
#my $params=$p->parse_param();	
#
#my $database=$params->{'database'};
#my $decoy_database=$database . "_decoy";
#
#if($database=~/,/)
#{
#	$decoy_database=~s/,/_decoy,/;
#}
#if($params->{'mass_correction}'}==1)
#{
#	my $mass_shift = retrieve("/home/jcho/dev/JUMPm/test/test_v1.8.0_20191030/4Plex_HILIC_Neg_3/4Plex_HILIC_Neg_3.1/.mass_shift");
#}
#
#my $falsified_shift = $params->{'mass_shift'};
#
#my @dtafiles = @ARGV;
#my $decoy_H_num = 1;
#if(defined($params->{'decoy_strategy'}))
#{
#	$decoy_H_num = $params->{'decoy_strategy'} - 1;
#}
#my $ms_hash_mol;



my $NC = 0;
my $CC = 0;
my $NC_std = 1;
my $CC_std = 1;
my $NC_defect_loc = 1.00335;
my $CC_defect_loc = 0.99703;
my $NC_defect_scale = 1;
my $CC_defect_scale = 1;



my $parameter = "jumpm_positive.params";
#my @dtafiles = ("/Research/Projects/7Metabolomics/JUMPm_v1.8.0_Xwang/Yeast_4Plex_RP_Pos_1/Yeast_4Plex_RP_Pos_1.1/14900.MS1");
#my @dtafiles = ("/home/xwang4/JUMPm/Yeast_4Plex_RP_Pos_1/Yeast_4Plex_RP_Pos_1.18/13002.MS1");
#my @dtafiles = glob("/Research/Projects/7Metabolomics/JUMPm_v1.8.0_Xwang/yeast_RP_POS_ph3_HM_01/yeast_RP_POS_ph3_HM_01.1/3*.MS1");
#my @dtafiles = ("/home/bxie/yeast4plex/yeast_RP_POS_ph3_HM_01/yeast_RP_POS_ph3_HM_01.1/30015.MS1");
my @dtafiles = ("/home/bxie/yeast4plex/yeast_RP_POS_ph3_HM_01/yeast_RP_POS_ph3_HM_01.1/1014.MS1");
my $p = Spiders::Params->new('-path'=>$parameter);
my $params = $p->parse_param();
my $database = $params->{'database'};
my $decoy_database = $database . "_decoy";
if ($database =~ /,/) {
	$decoy_database =~ s/,/_decoy,/;
}
if($params->{'mass_correction}'} == 1) {
	my $mass_shift = retrieve("/home/jcho/dev/JUMPm/test/test_v1.8.0_20191030/4Plex_HILIC_Neg_3/4Plex_HILIC_Neg_3.1/.mass_shift");
}
my $falsified_shift = $params->{'mass_shift'};
my $decoy_H_num = 1;
if(defined($params->{'decoy_strategy'})) {
	$decoy_H_num = $params->{'decoy_strategy'} - 1;
}
my $ms_hash_mol;





foreach(@dtafiles)
{
    my $index=1;
    my $dtafile = abs_path($_);
	print "\nProcessing scan: $dtafile\n";	
	
    my $scan = $dtafile;
    $scan =~s/(.*)\.(\d+)\.dta/$2/;	
    open(RESULT,">${dtafile}.smout");
    open(MRESULT,">${dtafile}.missile");	
	print MRESULT "Index","\t","C12 mass","\t","N15 mass","\t","C13 mass","\t","pair_score","\t","C12_C13_diff_score","\t", "C12_C13_int_rel","\t", "C12_C13_int_sim","\t","C12_N15_diff","\t", "C12_N15_int_rel","\t", "C12_N15_int_sim", "\t","Scan number", "\t","C12 intensity","\t","N15 intensity","\t","C13 intensity","\tType\tFormula\tAdduct\n";
		
   # print RESULT "MS1 scan: $dtafile\n";
    print RESULT "Index\tC12 mass\tN15 mass\tC13 mass\tMS1 Scan\tC12 Intensity\tN15 Intensity\tC13 Intensity\tPscore\tFormula\tAdduct\tName\tStructure (SMILES)\tInChI\tType\n";
	
	print "MS1 Deisotoping\n";			
	my $deisotope = new Spiders::Deisotope();  
	$deisotope->set_parameter($params);	
                        
    $deisotope->set_dta($dtafile);
	my ($mz_hash_ref,$int_ratio_hash,$peak_C,$charge_hash) = $deisotope->MS1_deisotope();

    #my $dta=new Spiders::Dta;
    #$dta->set_dta_file($dtafile) ;
    #$dta->process_dtafile($dtafile);

    #my %mz_hash = %{$dta->{'_mz_int_hash'}};
	%mz_hash = %$mz_hash_ref;	
    print scalar(keys %mz_hash)," features were detected\n";		
	
	print "Peaks pairing\n";			
    my $cluster;
	my $pair = new Spiders::Pairing();
	$pair->set_parameter($params);
	$pair->set_library_path("/home/jcho/dev/JUMPm/JUMPm_v1.8.0/R-3.1.0/");
	$pair->set_NC_ratio($NC);
	$pair->set_CC_ratio($CC);
	$pair->set_NC_std($NC_std);
	$pair->set_CC_std($CC_std);
	$pair->set_NC_defect_loc($NC_defect_loc);
	$pair->set_CC_defect_loc($CC_defect_loc);
	$pair->set_NC_defect_scale($NC_defect_scale);
	$pair->set_CC_defect_scale($CC_defect_scale);
	
    foreach my $mz (reverse sort {$mz_hash{$a}<=>$mz_hash{$b}} keys %mz_hash)
    {
        next if (!defined($mz_hash{$mz}));
		$cluster->{$mz}->{0}->{'mz'} = $mz;
        $cluster->{$mz}->{0}->{'int'} = $mz_hash{$mz};

        $pair->clustering($cluster,\%mz_hash,$mz);
    }
	
=head	
############### only for summary purpose ########
    my $dta=new Spiders::Dta;                  #
	#$dtafile .= ".iso";                        #
	$dta->set_prec_mz("1");                    #
    $dta->set_charge("0");	                    #
    $dta->set_dta_file($dtafile) ;	        #
	$dta->set_mz_int_hash(\%mz_hash);	
	$dta->write_dta_file();                    #
#################################################
=cut	
	
	
	my $dirname  = dirname($dtafile);
	
    my $cluster_mono = $pair->select_mono_peaks($cluster);

    my ($formula_comb) = $pair->pairing($cluster_mono,$dirname,$dtafile,$peak_C,$charge_hash);
	
	$index=1;	
	
	print scalar(keys %$formula_comb)," MISSILES were detected\n";
    foreach my $mono_mass (sort {$a<=>$b} keys %$formula_comb)
    {
	
		my $int_ratio = $int_ratio_hash->{$mono_mass};
        print "
performing mass formula database search for mass: $mono_mass\n";
		
        my @missile_formula = @{$formula_comb->{$mono_mass}->{'number'}};
        next if($mono_mass>1500);
		next if($#missile_formula<0);
###### For decoy testing		
		my $shift = 0;
		if($params->{'mode'} == -1)
		{
			$mono_mass = $mono_mass + 1.007277 + $shift;
		}
		else
		{
			$mono_mass = $mono_mass - 1.007277 + $shift;
		}
		if($params->{mass_correction})
		{
			$mono_mass = $mono_mass - $mono_mass * $mass_shift->{$scan} / 1000000;
		}		
		$mono_mass = $mono_mass + $falsified_shift;
		my $query = new Spiders::DatabaseQuery();
		$query->setBin("/home/jcho/dev/JUMPm/JUMPm_v1.8.0");		
		my ($return,$return_norm,$return_decoy) = $query->QueryMassWithDatabaseFilter($mono_mass, $params->{formula_mass_tolerance_searching}, $params->{mass_formula_database}, "yes",$database,$decoy_database);		
		# my $return =$query->QueryMassWithDatabase($mono_mass, $params->{formula_mass_tolerance}, $params->{mass_formula_database}, "yes");
		# my $return_decoy =$query->QueryMassWithDatabase($mono_mass, $params->{formula_mass_tolerance}, $params->{decoy_str_mass_formula_database}, "all");
		# my $return_norm =$query->QueryMassWithDatabase($mono_mass, $params->{formula_mass_tolerance}, $params->{norm_str_mass_formula_database}, "all");
		
		if($params->{'mode'} == -1)
		{
			$mono_mass = $mono_mass - 1.007277;
		}
		else
		{
			$mono_mass = $mono_mass + 1.007277;
		}
		
##### if the case of "mass only" 
		if($params->{'labeled_ID_method'} == 1)
		{
			foreach my $theoretical_formula (@$return) 
			{				
				chomp $theoretical_formula;
				my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);

				print MRESULT $index,"\t",$mono_mass,"\t","0","\t","0","\t","0","\t","0","\t","0","\t", "0","\t", "0","\t","0","\t", "0","\t", $scan, "\t", $formula_comb->{$mono_mass}->{'intensity'}, "\t","0","\t","0","\t", "pass", "\t", $formula_db, "\t \n";					
				$index++;					
			}
									
			foreach my $theoretical_formula (@$return_norm) 
			{				
				chomp $theoretical_formula;
				my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);
				print MRESULT $index,"\t",$mono_mass,"\t","0","\t","0","\t","0","\t","0","\t","0","\t", "0","\t", "0","\t","0","\t", "0","\t", $scan, "\t", $formula_comb->{$mono_mass}->{'intensity'}, "\t","0","\t","0","\t", "Target", "\t", $formula_db, "\t \n";

									
				$query->setDatabase($params->{structure_database});
				$query->setFormula($formula_db);
				$query->setMass($mass_db);
				#$query->setDatatype("ID,INCHIKEY,SMILES,GENERALNAME");

				$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
				#$query->setDatatype("SMILES");
				$query->setDBname($params->{database});				
				$str_db_return = $query->QueryStructureDatabase();				
	
				if(scalar(@$str_db_return)>0)
				{
					foreach (@$str_db_return)
					{
						chomp $_;					
						my ($ID,$formula_db,$InChI,$smiles,$IUPAC)=split(/\t/,$_);
						print RESULT $index,"\t",$mono_mass,"\t","0","\t","0","\t",$scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t","0","\t","0","\t","0","\t",$formula_db,"\t","0","\t",$IUPAC,"\t",$smiles,"\t",$InChI,"\tTarget","\n";						
					}
				}
						
				$index++;					
			}


			foreach my $theoretical_formula (@$return_decoy) 
			{				
				chomp $theoretical_formula;
				my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);
	

				if($formula_db=~/(.*H)(\d+)(\D+.*)/)
				{
					my $Hminus = $2 - 1  * $decoy_H_num;
					$updated_formula_db = $1 . $Hminus . $3
				}
	
							
				my $updated_mass_db = $mass_db - 1.00782503207  * $decoy_H_num;				
								
				my $query = new Spiders::DatabaseQuery();
				$query->setBin("/home/jcho/dev/JUMPm/JUMPm_v1.8.0");
				$query->setDatabase($params->{structure_database});
				$query->setFormula($updated_formula_db);
				$query->setMass($updated_mass_db);
				$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
				#$query->setDatatype("SMILES");						
				$query->setDBname($params->{database});
				
				my $str_db_return = $query->QueryStructureDatabase();	
				my $adduct = ($adduct_formula > 0) ? $adduct_name : "";
				$formula_db=~s/(\D+)1(\D+)/$1$2/g;
				$formula_db=~s/(\D+)1$/$1/g;									
				my $adduct_mass = ($adduct_formula > 0) ? $adduct_hash->{$adduct_name} : 0;                    						
				print MRESULT $index,"\t",$mono_mass,"\t",0,"\t",0,"\t",0,"\t",0,"\t",0,"\t", 0,"\t", 0,"\t",0,"\t", 0,"\t", $scan, "\t", $formula_comb->{$mono_mass}->{'intensity'}, "\t",0,"\t",0,"\t", "Decoy", "\t", $formula_db, "\t \n";	
								
				if(scalar(@$str_db_return)>0)
				{
					foreach (@$str_db_return)
					{
						chomp $_;
						my ($formula_db,$InChI,$smiles,$IUPAC)=split(/\t/,$_);

						print RESULT $index,"\t",$mono_mass,"\t","0","\t","0","\t",$scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t","0","\t","0","\t","0","\t",$formula_db,"\t", "0","\t",$IUPAC,"\t",$smiles,"\t",$InChI,"\tDecoy","\n";

						$index++;

					}
				}									
			}			
						
		}
		else
		{
			my $detected_formula=0;
			my $adduct_formula=0;			
			foreach my $missile_formula (@missile_formula)
			{
				my @N_C = split(/\,/,$missile_formula);
				my ($C_num_mz,$c_diff_score,$c_int_score, $c_similarity) = split(/\|/,$N_C[0]);
				my ($N_num_mz,$n_diff_score,$n_int_score, $n_similarity) = split(/\|/,$N_C[1]);
				my $pair_score = $N_C[2];
				my ($C_num,$C_mz) = split(/\:/,$C_num_mz);
				my ($N_num,$N_mz) = split(/\:/,$N_num_mz);	

###### correct mass before the search 		
				if($params->{mass_correction})
				{
					$C_mz = $C_mz - $C_mz * $mass_shift->{$scan} / 1000000;
					$N_mz = $N_mz - $N_mz * $mass_shift->{$scan} / 1000000;					
				}
				
					
				foreach my $theoretical_formula (@$return) 
				{
					chomp $theoretical_formula;
					my $carbon_match = 0;
					my $nitrogen_match = 0;
					
					if($params->{'labeled_ID_method'} == 2)
					{
						$carbon_match = ($theoretical_formula=~/$C_num\D+/);
						$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
						if($N_num eq "N0")
						{
							if($theoretical_formula!~/N/)
							{
								$nitrogen_match=1;
							}						
						}
					}
					elsif($params->{'labeled_ID_method'} == 3)
					{
						$carbon_match = ($theoretical_formula=~/$C_num\D+/);
						$nitrogen_match = 1;
					}
					elsif($params->{'labeled_ID_method'} == 4)
					{
						$carbon_match = 1;
						$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\
";
						exit(1);
					}					
										
					if($carbon_match == 1 and $nitrogen_match ==1)
					{
						
						my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);
						$formula_db=~s/(\D+)1(\D+)/$1$2/g;
						$formula_db=~s/(\D+)1$/$1/g;						
						
						print MRESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$pair_score,"\t",$c_diff_score,"\t",$c_int_score,"\t", $c_similarity,"\t", $n_diff_score,"\t",$n_int_score,"\t", $n_similarity,"\t", $scan, "\t", $formula_comb->{$mono_mass}->{'intensity'}, "\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t", "pass", "\t", $formula_db, "\t \n";
						$detected_formula++;						
					}
								
				}
				if($detected_formula==0 and $params->{'adduct'} == 1)
				{
					my $adduct_hash = $p->get_adducts($params);
					foreach $adduct_name (keys %$adduct_hash)
					{
						my $mono_mass_no_adduct = $mono_mass - $adduct_hash->{$adduct_name};
						my $N_mz_no_adduct = $C_mz_no_adduct = 0;
						if($N_mz>0)
						{
							$N_mz_no_adduct = $N_mz - $adduct_hash->{$adduct_name};
						}
						if($C_mz>0)
						{
							$C_mz_no_adduct = $C_mz - $adduct_hash->{$adduct_name};
						}					
						if($params->{'mode'} == -1)
						{
							$mono_mass_no_adduct = $mono_mass_no_adduct + 1.007277;
						}
						else
						{
							$mono_mass_no_adduct = $mono_mass_no_adduct - 1.007277;
						}					
						($return,$return_norm,$return_decoy) = $query->QueryMassWithDatabaseFilter($mono_mass_no_adduct, $params->{formula_mass_tolerance_searching}, $params->{mass_formula_database}, "yes",$database,$decoy_database);	
						if($params->{'mode'} == -1)
						{
							$mono_mass_no_adduct = $mono_mass_no_adduct - 1.007277;
						}
						else
						{
							$mono_mass_no_adduct = $mono_mass_no_adduct + 1.007277;
						}						
						
						foreach my $theoretical_formula (@$return) 
						{
							chomp $theoretical_formula;
							my $carbon_match = 0;
							my $nitrogen_match = 0;
							
							if($params->{'labeled_ID_method'} == 2)
							{
								$carbon_match = ($theoretical_formula=~/$C_num\D+/);
								$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
								if($N_num eq "N0")
								{
									if($theoretical_formula!~/N/)
									{
										$nitrogen_match=1;
									}						
								}
							}
							elsif($params->{'labeled_ID_method'} == 3)
							{
								$carbon_match = ($theoretical_formula=~/$C_num\D+/);
								$nitrogen_match = 1;
							}
							elsif($params->{'labeled_ID_method'} == 4)
							{
								$carbon_match = 1;
								$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\
";
								exit(1);
							}					
												
							if($carbon_match == 1 and $nitrogen_match ==1)
							{
								
								my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);
								my $adduct_mass = $adduct_hash->{$adduct_name};
								print MRESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$pair_score,"\t",$c_diff_score,"\t",$c_int_score,"\t", $c_similarity,"\t", $n_diff_score,"\t",$n_int_score,"\t", $n_similarity,"\t", $scan, "\t", $formula_comb->{$mono_mass}->{'intensity'}, "\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t", "pass", "\t", $formula_db, "\t",  $adduct_name,"\n";
								$adduct_formula++;						
							}
										
						}
						
						foreach my $theoretical_formula (@$return_norm) 
						{
						
						
							chomp $theoretical_formula;				
							my $carbon_match = 0;
							my $nitrogen_match = 0;
							
							if($params->{'labeled_ID_method'} == 2)
							{
								$carbon_match = ($theoretical_formula=~/$C_num\D+/);
								$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
								if($N_num eq "N0")
								{
									if($theoretical_formula!~/N/)
									{
										$nitrogen_match=1;
									}						
								}
							}
							elsif($params->{'labeled_ID_method'} == 3)
							{
								$carbon_match = ($theoretical_formula=~/$C_num\D+/);
								$nitrogen_match = 1;
							}
							elsif($params->{'labeled_ID_method'} == 4)
							{
								$carbon_match = 1;
								$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\n";
								exit(1);
							}					
								
							if($carbon_match == 1 and $nitrogen_match ==1)
							{
							
								my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);

								my $query = new Spiders::DatabaseQuery();
								$query->setBin("/home/jcho/dev/JUMPm/JUMPm_v1.8.0");										
								$query->setDatabase($params->{structure_database});
								$query->setFormula($formula_db);
								$query->setMass($mass_db);
								$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
								#$query->setDatatype("SMILES");
								$query->setDBname($params->{database});						
								$str_db_return = $query->QueryStructureDatabase();				
								my $adduct = ($adduct_formula > 0) ? $adduct_name : "";
								my $adduct_mass = ($adduct_formula > 0) ? $adduct_hash->{$adduct_name} : 0;
								$formula_db=~s/(\D+)1(\D+)/$1$2/g;
								$formula_db=~s/(\D+)1$/$1/g;									
								print MRESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$pair_score,"\t",$c_diff_score,"\t",$c_int_score,"\t", $c_similarity,"\t", $n_diff_score,"\t",$n_int_score,"\t", $n_similarity,"\t", $scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t","Target","\t",$formula_db, "\t",  $adduct,"\n";					
								
								if(scalar(@$str_db_return)>0)
								{
									foreach (@$str_db_return)
									{
										chomp $_;							
										my ($formula_db,$InChI,$smiles,$IUPAC)=split(/\t/,$_);

										print RESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t",$pair_score,"\t",$formula_db,"\t",  $adduct,"\t",$IUPAC,"\t",$smiles,"\t",$InChI,"\tTarget","\n";								
										
									}
								}
								
								$index++;		
							}	
						}
			################################################

						
						
			###### generate decoys 
						foreach my $theoretical_formula (@$return_decoy) 
						{
							chomp $theoretical_formula;
							my $carbon_match = 0;
							my $nitrogen_match = 0;
							
							if($params->{'labeled_ID_method'} == 2)
							{
								$carbon_match = ($theoretical_formula=~/$C_num\D+/);
								$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
								if($N_num eq "N0")
								{
									if($theoretical_formula!~/N/)
									{
										$nitrogen_match=1;
									}						
								}
							}
							elsif($params->{'labeled_ID_method'} == 3)
							{
								$carbon_match = ($theoretical_formula=~/$C_num\D+/);
								$nitrogen_match = 1;
							}
							elsif($params->{'labeled_ID_method'} == 4)
							{
								$carbon_match = 1;
								$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\n";
								exit(1);
							}					
								
							if($carbon_match == 1 and $nitrogen_match ==1)
							{				
								
								my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);

								my $updated_formula_db = $formula_db;
								if($formula_db=~/(.*H)(\d+)(\D+.*)/)
								{
									my $Hminus = $2 - 1  * $decoy_H_num;
									$updated_formula_db = $1 . $Hminus . $3
								}

											
								my $updated_mass_db = $mass_db - 1.00782503207  * $decoy_H_num;
									
								my $query = new Spiders::DatabaseQuery();
								$query->setBin("/home/jcho/dev/JUMPm/JUMPm_v1.8.0");
								$query->setDatabase($params->{structure_database});
								$query->setFormula($updated_formula_db);
								$query->setMass($updated_mass_db);
								$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
								#$query->setDatatype("SMILES");						
								$query->setDBname($params->{database});
								
								my $str_db_return = $query->QueryStructureDatabase();	
								my $adduct = ($adduct_formula > 0) ? $adduct_name : "";
								$formula_db=~s/(\D+)1(\D+)/$1$2/g;
								$formula_db=~s/(\D+)1$/$1/g;									
								my $adduct_mass = ($adduct_formula > 0) ? $adduct_hash->{$adduct_name} : 0;                    						
								print MRESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$pair_score,"\t",$c_diff_score,"\t",$c_int_score,"\t", $c_similarity,"\t", $n_diff_score,"\t",$n_int_score,"\t", $n_similarity,"\t", $scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t","Decoy","\t",$formula_db, "\t",  $adduct,"\n";
								
								if(scalar(@$str_db_return)>0)
								{
									foreach (@$str_db_return)
									{
										chomp $_;
										my ($formula_db,$InChI,$smiles,$IUPAC)=split(/\t/,$_);

										print RESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t",$pair_score,"\t",$formula_db,"\t",  $adduct,"\t",$IUPAC,"\t",$smiles,"\t",$InChI,"\tDecoy","\n";

										my $smiles_orig = $smiles;
										$index++;

									}
								}						
								
								$index++;		
							}	
						}
					}
	
				}
				else
				{

####### generating structure and scoring #######################
			
					foreach my $theoretical_formula (@$return_norm) 
					{
					
					
						chomp $theoretical_formula;				
						my $carbon_match = 0;
						my $nitrogen_match = 0;
						
						if($params->{'labeled_ID_method'} == 2)
						{
							$carbon_match = ($theoretical_formula=~/$C_num\D+/);
							$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
							if($N_num eq "N0")
							{
								if($theoretical_formula!~/N/)
								{
									$nitrogen_match=1;
								}						
							}
						}
						elsif($params->{'labeled_ID_method'} == 3)
						{
							$carbon_match = ($theoretical_formula=~/$C_num\D+/);
							$nitrogen_match = 1;
						}
						elsif($params->{'labeled_ID_method'} == 4)
						{
							$carbon_match = 1;
							$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
						}
						else
						{
							print "please set the labeled_ID_method parameter [1-4]\n";
							exit(1);
						}					
							
						if($carbon_match == 1 and $nitrogen_match ==1)
						{
						
							my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);

							my $query = new Spiders::DatabaseQuery();
							$query->setBin("/home/jcho/dev/JUMPm/JUMPm_v1.8.0");										
							$query->setDatabase($params->{structure_database});
							$query->setFormula($formula_db);
							$query->setMass($mass_db);
							$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
							#$query->setDatatype("SMILES");
							$query->setDBname($params->{database});						
							$str_db_return = $query->QueryStructureDatabase();				

							$formula_db=~s/(\D+)1(\D+)/$1$2/g;
							$formula_db=~s/(\D+)1$/$1/g;								
							print MRESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$pair_score,"\t",$c_diff_score,"\t",$c_int_score,"\t", $c_similarity,"\t", $n_diff_score,"\t",$n_int_score,"\t", $n_similarity,"\t", $scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t","Target","\t",$formula_db, "\t",  $adduct,"\n";					
							
							if(scalar(@$str_db_return)>0)
							{
								foreach (@$str_db_return)
								{
									chomp $_;							
									my ($formula_db,$InChI,$smiles,$IUPAC)=split(/\t/,$_);

									print RESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t",$pair_score,"\t",$formula_db,"\t",  "","\t",$IUPAC,"\t",$smiles,"\t",$InChI,"\tTarget","\n";								
									
								}
							}
							
							$index++;		
						}	
					}
		################################################

					
					
		###### generate decoys 
					foreach my $theoretical_formula (@$return_decoy) 
					{
						chomp $theoretical_formula;
						my $carbon_match = 0;
						my $nitrogen_match = 0;
						
						if($params->{'labeled_ID_method'} == 2)
						{
							$carbon_match = ($theoretical_formula=~/$C_num\D+/);
							$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
							if($N_num eq "N0")
							{
								if($theoretical_formula!~/N/)
								{
									$nitrogen_match=1;
								}						
							}
						}
						elsif($params->{'labeled_ID_method'} == 3)
						{
							$carbon_match = ($theoretical_formula=~/$C_num\D+/);
							$nitrogen_match = 1;
						}
						elsif($params->{'labeled_ID_method'} == 4)
						{
							$carbon_match = 1;
							$nitrogen_match = ($theoretical_formula=~/$N_num\D+/);
						}
						else
						{
							print "please set the labeled_ID_method parameter [1-4]\n";
							exit(1);
						}					
							
						if($carbon_match == 1 and $nitrogen_match ==1)
						{				
							
							my ($formula_db,$mass_db) = split(/:/,$theoretical_formula);

							my $updated_formula_db = $formula_db;
							if($formula_db=~/(.*H)(\d+)(\D+.*)/)
							{
								my $Hminus = $2 - 1  * $decoy_H_num;
								$updated_formula_db = $1 . $Hminus . $3
							}

										
							my $updated_mass_db = $mass_db - 1.00782503207  * $decoy_H_num;

							my $query = new Spiders::DatabaseQuery();
							$query->setBin("/home/jcho/dev/JUMPm/JUMPm_v1.8.0");
							$query->setDatabase($params->{structure_database});
							$query->setFormula($updated_formula_db);
							$query->setMass($updated_mass_db);
							$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
							#$query->setDatatype("SMILES");						
							$query->setDBname($params->{database});
							
							my $str_db_return = $query->QueryStructureDatabase();	
							$formula_db=~s/(\D+)1(\D+)/$1$2/g;
							$formula_db=~s/(\D+)1$/$1/g;	                  						
							print MRESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$pair_score,"\t",$c_diff_score,"\t",$c_int_score,"\t", $c_similarity,"\t", $n_diff_score,"\t",$n_int_score,"\t", $n_similarity,"\t", $scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t","Decoy","\t",$formula_db, "\t",  $adduct,"\n";
							
							if(scalar(@$str_db_return)>0)
							{
								foreach (@$str_db_return)
								{
									chomp $_;
									my ($formula_db,$InChI,$smiles,$IUPAC)=split(/\t/,$_);

									print RESULT $index,"\t",$mono_mass,"\t",$N_mz,"\t",$C_mz,"\t",$scan,"\t",$formula_comb->{$mono_mass}->{'intensity'},"\t",$formula_comb->{$N_mz}->{'intensity'},"\t",$formula_comb->{$C_mz}->{'intensity'},"\t",$pair_score,"\t",$formula_db,"\t",  "","\t",$IUPAC,"\t",$smiles,"\t",$InChI,"\tDecoy","\n";

									my $smiles_orig = $smiles;
									$index++;

								}
							}						
							
							$index++;		
						}	
					}
				}
			}			
		}	
    }
	
}
		
