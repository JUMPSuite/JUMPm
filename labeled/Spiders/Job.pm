#!/usr/bin/perl

######### Job #################################################
#                                                             #
#       **************************************************    #  
#       **** Job module for JUMPm   		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2015 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::Job;


use warnings;
use File::Basename;
use Storable;
use Parallel::ForkManager;
use Cwd;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.03;


@ISA	 = qw(Exporter);
@EXPORT = qw(set_dta_path get_dta_path set_pip get_pip set_library_path get_library_path set_R_library_path get_R_library_path set_parameter get_parameter set_mass_shift get_mass_shift set_log_file get_log_file create_frag_script create_pscript create_mscript_partial carbon_max_number create_mscript_partial2 carbon_max_number create_mscript_unlabel carbon_number_matching create_mscript Check_Job_stat LuchParallelJob runjobs luanch_ms2_jobs);
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
        _dta_path => undef,
        _sim_path  => undef,
    };
    bless $self, $class;
    return $self;
}

sub set_dta_path
{
	my ($self,$dta_path)=@_;
	$self->{_dta_path}=$dta_path;
}

sub get_dta_path
{
	my $self=shift;
	return $self->{_dta_path};
}

sub set_pip
{
	my ($self,$pip)=@_;
	$self->{_pip}=$pip;	
}

sub get_pip
{
	my $self=shift;
	return $self->{_pip};
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

sub set_R_library_path
{
	my ($self,$lib)=@_;
	$self->{_R_lib}=$lib;	
}

sub get_R_library_path
{
	my ($self)=@_;
	return $self->{_R_lib};	
}


sub set_parameter
{
	my ($self,$param)=@_;
	$self->{'_parameter'}=$param;
}


sub get_parameter
{
	my $self=shift;
	return $self->{'_parameter'};
}

sub set_mass_shift
{
	my ($self,$mass_shift)=@_;
	$self->{'_mass_shift'}=$mass_shift;
}


sub get_mass_shift
{
	my $self=shift;
	return $self->{'_mass_shift'};
}

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

sub create_frag_script
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	my $parameter = $self->get_parameter();
	
	my $neutralLossFile = $lib . "/neutralLoss.csv";	
	my $bondEnergyFile =  $lib . "/bondenergies.txt";


	
	open(RUNSHELL,">$dir/frag_shell.pl");
=head
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::MS2_scoring;
use Spiders::Dta;

my ($help,$dtafile,$smile,$mass,$depth);
GetOptions('-help|h'=>\\\$help,
		'-dtafile=s'=>\\\$dtafile,
		'-smile=s'=>\\\$smile,
		'-mass=s'=>\\\$mass,
		'-depth=s'=>\\\$depth,	
		);
		
		SmileFragmentation(\$smile, \$mass, \$depth, \$neutralLossFile, \$bondEnergyFile);
		my \@frag_mz_values=();
		for(my \$i=2;\$i<\$\#frag_return;\$i++) 
		{
			my \@frag_data = split(/\\s+/,\$frag_return[\$i]);
			next if(\$frag_data[-1]=~/[CHNOPS]/);
			push (\@frag_mz_values,\$frag_data[-1]);
		}
		my \$dta=new Spiders::Dta;  
		\$dta->set_dta_file(\$dtafile) ;
		\$dta->process_dtafile(\$dtafile);

		my \%mz_hash = %{\$dta->{'_mz_hash'}};
	
		my \$scoring = new Spiders::MS2_scoring;
		\$scoring->compare_theoritical_experiment(\\\%mz_hash,\\\@frag_mz_values);
		\$scoring->write_result();
EOF
=cut
	close(RUNSHELL);			
}

sub create_pscript
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();

	
	open(PAIRSHELL,">$dir/pair_shell.pl");
print PAIRSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,		
		);
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	
my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

my \@dtafiles = \@ARGV;


my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;

    my \$dtafile = abs_path(\$_);


	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\)\\\.\(\\d\+\)\\\.dta\/\$2\/;	
    open(MRESULT,">\${dtafile}.pairs");	
    open(MRESULTALL,">\${dtafile}.pairs.pass");		
    my \$dta=new Spiders::Dta;  
    \$dta->set_dta_file(\$dtafile) ;
    \$dta->process_dtafile(\$dtafile);

    my \%mz_hash = %{\$dta->{'_mz_int_hash'}};

	
    my \$cluster;
	my \$pair = new Spiders::Pairing();
	\$pair->set_parameter(\$params);
	\$pair->set_NC_ratio(0);
	\$pair->set_CC_ratio(0);	
    foreach my \$mz (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
            next if (\!defined(\$mz_hash{\$mz}));
			\$cluster->{\$mz}->{0}->{'mz'} = \$mz;
            \$cluster->{\$mz}->{0}->{'int'} = \$mz_hash{\$mz};

            \$pair->clustering(\$cluster,\\\%mz_hash,\$mz);
    }
	
	my \$dirname  = dirname(\$dtafile);
    my \$cluster_mono = \$pair->select_mono_peaks(\$cluster);
    my (\$formula_comb) = \$pair->init_pairing(\$cluster_mono,\$dirname);
	
    foreach my \$mono_mass (sort {\$a<=>\$b} keys \%\$formula_comb)
    {
        my \@missile_formula = \@{\$formula_comb->{\$mono_mass}->{'number'}};
        next if(\$mono_mass>1500);
		next if(\$#missile_formula<0);
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass + 1.007277 + \$shift;
		}
		else
		{
			\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		}
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");		
		my (\$return_all,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_pairing}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);
		#my \$return_norm = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all", "hmdb");
		#my \$return_norm = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all", "pubchem");
		#my \$return_norm = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all", "kegg");
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass - 1.007277;
		}
		else
		{
			\$mono_mass = \$mono_mass + 1.007277;
		}

		foreach my \$missile_formula (\@missile_formula)
		{
			my \@N_C = split(/\\,/,\$missile_formula);
			my (\$C_num_mz,\$C_mass_defect,\$C_relative_int) = split(/\\|/,\$N_C[0]);
			my (\$N_num_mz,\$N_mass_defect,\$N_relative_int) = split(/\\|/,\$N_C[1]);
				
			my (\$C_num,\$C_mz) = split(/\\:/,\$C_num_mz);
			my (\$N_num,\$N_mz) = split(/\\:/,\$N_num_mz);	
					
			foreach my \$theoretical_formula (\@\$return_norm) 
			{
				chomp \$theoretical_formula;
				my \$carbon_match = 0;
				my \$nitrogen_match = 0;

				\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
				\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
=head				
				if(\$N_num eq "N0")
				{
					if(\$theoretical_formula!~/N/)
					{
						\$nitrogen_match=1;
					}						
				}				
=cut				
				if(\$carbon_match == 1 and \$nitrogen_match ==1)
				{				
					print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$C_mass_defect,"\\t", \$C_relative_int,"\\t", \$N_mass_defect,"\\t",\$N_relative_int,"\\n";					
				}
			}
		
			foreach my \$theoretical_formula (\@\$return_all) 
			{
				chomp \$theoretical_formula;
				my \$carbon_match = 0;
				my \$nitrogen_match = 0;

				\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
				\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
=head					
				if(\$N_num eq "N0")
				{
					if(\$theoretical_formula!~/N/)
					{
						\$nitrogen_match=1;
					}						
				}				
=cut					
				if(\$carbon_match == 1 and \$nitrogen_match ==1)
				{				
					print MRESULTALL \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$C_mass_defect,"\\t", \$C_relative_int,"\\t", \$N_mass_defect,"\\t",\$N_relative_int,"\\n";					
				}
			}
		
		}					
	}
}	
EOF
	close(PAIRSHELL);	
}


sub create_mscript_partial
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	my $R_lib = $self->get_R_library_path();

	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,
		'-NC=f'=>\\\$NC,
		'-CC=f'=>\\\$CC,
		'-NC_std=f'=>\\\$NC_std,
		'-CC_std=f'=>\\\$CC_std,
		'-NC_defect_loc=f'=>\\\$NC_defect_loc,
		'-CC_defect_loc=f'=>\\\$CC_defect_loc,
		'-NC_defect_scale=f'=>\\\$NC_defect_scale,
		'-CC_defect_scale=f'=>\\\$CC_defect_scale,
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

if(\$database=~/\,/)
{
	\$decoy_database=~s/\,/\_decoy\,/;
}
if(\$params->{'mass_correction}'}==1)
{
	my \$mass_shift = retrieve("$dir/.mass_shift");
}
my \@dtafiles = \@ARGV;
my \$decoy_H_num = 1;
if(defined(\$params->{'decoy_strategy'}))
{
	\$decoy_H_num = \$params->{'decoy_strategy'} - 1;
}
my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;
    my \$dtafile = abs_path(\$_);
	print "\\nProcessing scan: \$dtafile\\n";	
	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\)\\\.\(\\d\+\)\\\.dta\/\$2\/;	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_diff_score","\\t", "C12_C13_int_rel","\\t", "C12_C13_int_sim","\\t","C12_N15_diff","\\t", "C12_N15_int_rel","\\t", "C12_N15_int_sim", "\\t","Scan number", "\\t","C12 intensity","\\t","N15 intensity","\\t","C13 intensity","\\tType\\tFormula\\n";
		
   # print RESULT "MS1 scan: \$dtafile\\n";
    print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tPscore\\tFormula\\tName\\tStructure (SMILES)\\tInChI\\tType\\n";
	
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	my (\$mz_hash_ref,\$int_ratio_hash,\$peak_C,\$charge_hash) = \$deisotope->MS1_deisotope();
    #my \$dta=new Spiders::Dta;  
    #\$dta->set_dta_file(\$dtafile) ;
    #\$dta->process_dtafile(\$dtafile);

    #my \%mz_hash = %{\$dta->{'_mz_int_hash'}};
	\%mz_hash = \%\$mz_hash_ref;	
    print scalar(keys \%mz_hash)," features were detected\\n";		
	
	print "Peaks pairing\\n";			
    my \$cluster;
	my \$pair = new Spiders::Pairing();
	\$pair->set_parameter(\$params);
	\$pair->set_library_path("$R_lib");
	\$pair->set_NC_ratio(\$NC);
	\$pair->set_CC_ratio(\$CC);
	\$pair->set_NC_std(\$NC_std);
	\$pair->set_CC_std(\$CC_std);
	\$pair->set_NC_defect_loc(\$NC_defect_loc);
	\$pair->set_CC_defect_loc(\$CC_defect_loc);
	\$pair->set_NC_defect_scale(\$NC_defect_scale);
	\$pair->set_CC_defect_scale(\$CC_defect_scale);
	
    foreach my \$mz (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
        next if (\!defined(\$mz_hash{\$mz}));
		\$cluster->{\$mz}->{0}->{'mz'} = \$mz;
        \$cluster->{\$mz}->{0}->{'int'} = \$mz_hash{\$mz};

        \$pair->clustering(\$cluster,\\\%mz_hash,\$mz);
    }
	
	
############### only for summary purpose ########
    my \$dta=new Spiders::Dta;                  #
	#\$dtafile .= ".iso";                        #
	\$dta->set_prec_mz("1");                    #
    \$dta->set_charge("0");	                    #
    \$dta->set_dta_file(\$dtafile) ;	        #
	\$dta->set_mz_int_hash(\\\%mz_hash);	
	\$dta->write_dta_file();                    #
#################################################
	
	
	
	my \$dirname  = dirname(\$dtafile);
	
    my \$cluster_mono = \$pair->select_mono_peaks(\$cluster);

    my (\$formula_comb) = \$pair->pairing(\$cluster_mono,\$dirname,\$dtafile,\$peak_C,\$charge_hash);
	
	\$index=1;	
	
	print scalar(keys \%\$formula_comb)," MISSILES were detected\\n";
    foreach my \$mono_mass (sort {\$a<=>\$b} keys \%\$formula_comb)
    {
		my \$int_ratio = \$int_ratio_hash->{\$mono_mass};
        print "\nperforming mass formula database search for mass: \$mono_mass\\n";
		
        my \@missile_formula = \@{\$formula_comb->{\$mono_mass}->{'number'}};
        next if(\$mono_mass>1500);
		next if(\$#missile_formula<0);
###### For decoy testing		
		my \$shift = 0;
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass + 1.007277 + \$shift;
		}
		else
		{
			\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		}
		if(\$params->{mass_correction})
		{
			\$mono_mass = \$mono_mass - \$mono_mass * \$mass_shift->{\$scan} / 1000000;
		}		
	
		
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");		
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);		
		# my \$return =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{mass_formula_database}, "yes");
		# my \$return_decoy =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{decoy_str_mass_formula_database}, "all");
		# my \$return_norm =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all");
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass - 1.007277;
		}
		else
		{
			\$mono_mass = \$mono_mass + 1.007277;
		}
		
##### if the case of "mass only" 
		if(\$params->{'labeled_ID_method'} == 1)
		{
			foreach my \$theoretical_formula (\@\$return) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				
				\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
				\$formula_db=~s/(\\D+)1\$/\$1/g;					
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\n";
						
				\$index++;					
			}
									
			foreach my \$theoretical_formula (\@\$return_norm) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
				\$formula_db=~s/(\\D+)1\$/\$1/g;					
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\tnorm\\n";

									
				\$query->setDatabase(\$params->{structure_database});
				\$query->setFormula(\$formula_db);
				\$query->setMass(\$mass_db);
				\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
				#\$query->setDatatype("SMILES");
				\$query->setDBname(\$params->{database});				
				\$str_db_return = \$query->QueryStructureDatabase();				
	
				if(scalar(\@\$str_db_return)>0)
				{
					foreach (\@\$str_db_return)
					{
						chomp \$_;					
						my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
						print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$formula_db,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
					}
				}
						
				\$index++;					
			}


			foreach my \$theoretical_formula (\@\$return_decoy) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
				\$formula_db=~s/(\\D+)1\$/\$1/g;					
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\tdecoy\\n";
		
		
				\$index++;					
			}			
			
			
		}
		else
		{		
			foreach my \$missile_formula (\@missile_formula)
			{
				my \@N_C = split(/\\,/,\$missile_formula);
				my (\$C_num_mz,\$c_diff_score,\$c_int_score, \$c_similarity) = split(/\\|/,\$N_C[0]);
				my (\$N_num_mz,\$n_diff_score,\$n_int_score, \$n_similarity) = split(/\\|/,\$N_C[1]);
				my \$pair_score = \$N_C[2];
				my (\$C_num,\$C_mz) = split(/\\:/,\$C_num_mz);
				my (\$N_num,\$N_mz) = split(/\\:/,\$N_num_mz);	

###### correct mass before the search 		
				if(\$params->{mass_correction})
				{
					\$C_mz = \$C_mz - \$C_mz * \$mass_shift->{\$scan} / 1000000;
					\$N_mz = \$N_mz - \$N_mz * \$mass_shift->{\$scan} / 1000000;					
				}
				
					
				foreach my \$theoretical_formula (\@\$return) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
				
						my \$max_C_num = carbon_max_number(\$int_ratio);
											
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
												
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						#\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						my \$max_C_num = carbon_max_number(\$int_ratio);						
						\$nitrogen_match = 1;
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}							
						}						
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{
						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;							
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t", "pass", "\\t", \$formula_db, "\\n";			
					}
				}
			


####### generating structure and scoring #######################
			
				foreach my \$theoretical_formula (\@\$return_norm) 
				{
				
				
					chomp \$theoretical_formula;				
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
					
						my \$max_C_num = carbon_max_number(\$int_ratio);
											
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
												
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}
						}					
					
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						#\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						my \$max_C_num = carbon_max_number(\$int_ratio);						
						\$nitrogen_match = 1;
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}							
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{
					
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

						my \$query = new Spiders::DatabaseQuery();
						\$query->setBin("$lib");										
						\$query->setDatabase(\$params->{structure_database});
						\$query->setFormula(\$formula_db);
						\$query->setMass(\$mass_db);
						\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
						#\$query->setDatatype("SMILES");
						\$query->setDBname(\$params->{database});						
						\$str_db_return = \$query->QueryStructureDatabase();				
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;	
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Target","\\t",\$formula_db,"\\n";					
						
						if(scalar(\@\$str_db_return)>0)
						{
							foreach (\@\$str_db_return)
							{
								chomp \$_;							
								my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
								print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
							}
						}
						
						\$index++;		
					}	
				}
	################################################

				
				
	###### generate decoys 
				foreach my \$theoretical_formula (\@\$return_decoy) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
						my \$max_C_num = carbon_max_number(\$int_ratio);
											
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
												
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						#\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						my \$max_C_num = carbon_max_number(\$int_ratio);						
						\$nitrogen_match = 1;
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}							
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{				
						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

						my \$updated_formula_db = \$formula_db;
						if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
						{
							my \$Hminus = \$2 - \$decoy_H_num;
							\$updated_formula_db = \$1 . \$Hminus . \$3
						}

									
						my \$updated_mass_db = \$mass_db - 1.00782503207 * \$decoy_H_num;
					
						my \$query = new Spiders::DatabaseQuery();
						\$query->setBin("$lib");
						\$query->setDatabase(\$params->{structure_database});
						\$query->setFormula(\$updated_formula_db);
						\$query->setMass(\$updated_mass_db);
						\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
						#\$query->setDatatype("SMILES");						
						\$query->setDBname(\$params->{database});
						
						my \$str_db_return = \$query->QueryStructureDatabase();				
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;	
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Decoy","\\t",\$formula_db,"\\n";
						
						if(scalar(\@\$str_db_return)>0)
						{
							foreach (\@\$str_db_return)
							{
								chomp \$_;
								my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
								print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
								
								my \$smiles_orig = \$smiles;
								\$index++;

							}
						}						
						
						\$index++;		
					}	
				}
			}				
		}	
    }
	
}
	
sub carbon_max_number
{
	my (\$int_ratio) = \@_;
	my \$max_C_num = int(\$int_ratio + \$int_ratio*0.3)+1;		
	return \$max_C_num;
}	
	
EOF
	close(RUNSHELL);	
}


sub create_mscript_partial2
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	my $R_lib = $self->get_R_library_path();

	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,
		'-NC=f'=>\\\$NC,
		'-CC=f'=>\\\$CC,
		'-NC_std=f'=>\\\$NC_std,
		'-CC_std=f'=>\\\$CC_std,
		'-NC_defect_loc=f'=>\\\$NC_defect_loc,
		'-CC_defect_loc=f'=>\\\$CC_defect_loc,
		'-NC_defect_scale=f'=>\\\$NC_defect_scale,
		'-CC_defect_scale=f'=>\\\$CC_defect_scale,
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

if(\$database=~/\,/)
{
	\$decoy_database=~s/\,/\_decoy\,/;
}
if(\$params->{'mass_correction}'}==1)
{
	my \$mass_shift = retrieve("$dir/.mass_shift");
}

my \$falsified_shift = \$params->{'mass_shift'};

if(\$params->{'mass_correction}'}==1)
{
	my \$mass_shift = retrieve("$dir/.mass_shift");
}


my \@dtafiles = \@ARGV;
my \$decoy_H_num = 1;
if(defined(\$params->{'decoy_strategy'}))
{
	\$decoy_H_num = \$params->{'decoy_strategy'} - 1;
}
my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;
    my \$dtafile = abs_path(\$_);
	print "\\nProcessing scan: \$dtafile\\n";	
	
    my \$scan = basename(\$dtafile);
    \$scan =~s\/\(\\d\+\)\\\.MS1\/\$1\/;	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_diff_score","\\t", "C12_C13_int_rel","\\t", "C12_C13_int_sim","\\t","C12_N15_diff","\\t", "C12_N15_int_rel","\\t", "C12_N15_int_sim", "\\t","Scan number", "\\t","C12 intensity","\\t","N15 intensity","\\t","C13 intensity","\\tType\\tFormula\\n";
		
   # print RESULT "MS1 scan: \$dtafile\\n";
    print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tPscore\\tFormula\\tName\\tStructure (SMILES)\\tInChI\\tType\\n";
	
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	my (\$mz_hash_ref,\$int_ratio_hash,\$peak_C,\$charge_hash) = \$deisotope->MS1_deisotope();
    #my \$dta=new Spiders::Dta;  
    #\$dta->set_dta_file(\$dtafile) ;
    #\$dta->process_dtafile(\$dtafile);

    #my \%mz_hash = %{\$dta->{'_mz_int_hash'}};
	\%mz_hash = \%\$mz_hash_ref;	
    print scalar(keys \%mz_hash)," features were detected\\n";		
	
	print "Peaks pairing\\n";			
    my \$cluster;
	my \$pair = new Spiders::Pairing();
	\$pair->set_parameter(\$params);
	\$pair->set_library_path("$R_lib");
	\$pair->set_NC_ratio(\$NC);
	\$pair->set_CC_ratio(\$CC);
	\$pair->set_NC_std(\$NC_std);
	\$pair->set_CC_std(\$CC_std);
	\$pair->set_NC_defect_loc(\$NC_defect_loc);
	\$pair->set_CC_defect_loc(\$CC_defect_loc);
	\$pair->set_NC_defect_scale(\$NC_defect_scale);
	\$pair->set_CC_defect_scale(\$CC_defect_scale);
	
    foreach my \$mz (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
        next if (\!defined(\$mz_hash{\$mz}));
		\$cluster->{\$mz}->{0}->{'mz'} = \$mz;
        \$cluster->{\$mz}->{0}->{'int'} = \$mz_hash{\$mz};

        \$pair->clustering(\$cluster,\\\%mz_hash,\$mz);
    }
	
	
############### only for summary purpose ########
    my \$dta=new Spiders::Dta;                  #
	#\$dtafile .= ".iso";                        #
	\$dta->set_prec_mz("1");                    #
    \$dta->set_charge("0");	                    #
    \$dta->set_dta_file(\$dtafile) ;	        #
	\$dta->set_mz_int_hash(\\\%mz_hash);	
	\$dta->write_dta_file();                    #
#################################################
	
	
	
	my \$dirname  = dirname(\$dtafile);
	
    my \$cluster_mono = \$pair->select_mono_peaks(\$cluster);

    my (\$formula_comb) = \$pair->pairing(\$cluster_mono,\$dirname,\$dtafile,\$peak_C,\$charge_hash);
	
	\$index=1;	
	
	print scalar(keys \%\$formula_comb)," MISSILES were detected\\n";
    foreach my \$mono_mass (sort {\$a<=>\$b} keys \%\$formula_comb)
    {
		my \$int_ratio = \$int_ratio_hash->{\$mono_mass};
        print "\nperforming mass formula database search for mass: \$mono_mass\\n";
		
        my \@missile_formula = \@{\$formula_comb->{\$mono_mass}->{'number'}};
        next if(\$mono_mass>1500);
		next if(\$#missile_formula<0);
###### For decoy testing		
		my \$shift = 0;
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass + 1.007277 + \$shift;
		}
		else
		{
			\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		}
		if(\$params->{mass_correction})
		{
			\$mono_mass = \$mono_mass - \$mono_mass * \$mass_shift->{\$scan} / 1000000;
		}		
	
		\$mono_mass = \$mono_mass + \$falsified_shift;
		
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");
		\$query->set_parameter(\$params);
		#my \$mzXML_file = $dir;
		#\$mzXML_file =~s/(.*)\\.(\\d+)/\$1/;
		my \$basename = basename(\$dtafile);
		my \$dirname = dirname(\$dtafile);
		
		my \$dtafile_orignal = \$dirname . "\/MS1\/" . \$basename;
		
		
		\$query->set_dta_file(\"\$dtafile_orignal\");
		
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);		
		# my \$return =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{mass_formula_database}, "yes");
		# my \$return_decoy =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{decoy_str_mass_formula_database}, "all");
		# my \$return_norm =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all");
		
		
		
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass - 1.007277;
		}
		else
		{
			\$mono_mass = \$mono_mass + 1.007277;
		}
		
##### if the case of "mass only" 
		if(\$params->{'labeled_ID_method'} == 1)
		{
			foreach my \$theoretical_formula (\@\$return) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
			
				\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
				\$formula_db=~s/(\\D+)1\$/\$1/g;					
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\n";
						
				\$index++;					
			}
									
			foreach my \$theoretical_formula (\@\$return_norm) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
				\$formula_db=~s/(\\D+)1\$/\$1/g;					
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\tnorm\\n";

									
				\$query->setDatabase(\$params->{structure_database});
				\$query->setFormula(\$formula_db);
				\$query->setMass(\$mass_db);
				\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
				#\$query->setDatatype("SMILES");
				\$query->setDBname(\$params->{database});				
				\$str_db_return = \$query->QueryStructureDatabase();				
	
				if(scalar(\@\$str_db_return)>0)
				{
					foreach (\@\$str_db_return)
					{
						chomp \$_;					
						my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
						print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$formula_db,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
					}
				}
						
				\$index++;					
			}


			foreach my \$theoretical_formula (\@\$return_decoy) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\tdecoy\\n";
		
		
				\$index++;					
			}			
			
			
		}
		else
		{		
			foreach my \$missile_formula (\@missile_formula)
			{
				my \@N_C = split(/\\,/,\$missile_formula);
				my (\$C_num_mz,\$c_diff_score,\$c_int_score, \$c_similarity) = split(/\\|/,\$N_C[0]);
				my (\$N_num_mz,\$n_diff_score,\$n_int_score, \$n_similarity) = split(/\\|/,\$N_C[1]);
				my \$pair_score = \$N_C[2];
				my (\$C_num,\$C_mz) = split(/\\:/,\$C_num_mz);
				my (\$N_num,\$N_mz) = split(/\\:/,\$N_num_mz);	
				\$C_num =~s/C//;
###### correct mass before the search 		
				if(\$params->{mass_correction})
				{
					\$C_mz = \$C_mz - \$C_mz * \$mass_shift->{\$scan} / 1000000;
					\$N_mz = \$N_mz - \$N_mz * \$mass_shift->{\$scan} / 1000000;					
				}
				
				my \$formula_unlabel = "";	
				foreach my \$theoretical_formula (\@\$return) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
				
						my \$max_C_num = carbon_max_number(\$int_ratio);
											
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
												
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/C\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
									
						next if((\$C_num_mz/12)>5*\$C_num);
						next if(\$theoretical_formula=~/F/);
						next if(\$theoretical_formula=~/Br/);						
						#\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						
						my \$max_C_num = carbon_max_number(\$int_ratio);						
						\$nitrogen_match = 1;
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/C\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}							
						}						
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\\n";
						exit(1);
					}					
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{

						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
						\$formula_unlabel .= "\${formula_db},";
					}
					
				}
		
				if(length(\$formula_unlabel)>1)
				{
					my \$min_mass = \$C_mz - 5*1.007277;
					my \$max_mass = \$C_mz + 5*1.007277;				
					my \$value = \$query->partial_candidate(\$formula_unlabel, \$min_mass, \$max_mass);
					print \@\$value;
					foreach my \$result (\@\$value)
					{
						my \@data = split(/\\\t/,\$result);
						\$formula_db = \$data[0];
						my \$score = \$data[4];
						\$score=~s/MatchScore\://;						
						next if(\$score<0.9);
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;							
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t", "pass", "\\t", \$formula_db, "\\n";
					}
					
				}
####### generating structure and scoring #######################
				my \$formula_unlabel_normal = "";			
				foreach my \$theoretical_formula (\@\$return_norm) 
				{
				
				
					chomp \$theoretical_formula;				
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
					
						my \$max_C_num = carbon_max_number(\$int_ratio);
											
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
												
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/C\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}
						}					
					
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						#\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						my \$max_C_num = carbon_max_number(\$int_ratio);						
						\$nitrogen_match = 1;
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/C\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}							
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{
					
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
						\$formula_unlabel_normal .= "\${formula_db},";						
						
						


					}	
				}
				
				if(length(\$formula_unlabel_normal)>1)
				{
					my \$min_mass = \$C_mz - 5*1.007277;
					my \$max_mass = \$C_mz + 5*1.007277;				
					my \$value = \$query->partial_candidate(\$formula_unlabel_normal, \$min_mass, \$max_mass);
					print \@\$value;
					foreach my \$result (\@\$value)
					{
						my \@data = split(/\\\t/,\$result);
						\$formula_db = \$data[0];
						my \$score = \$data[4];
						\$score=~s/MatchScore\://;						
						next if(\$score<0.9);
						
						#print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t", "pass", "\\t", \$formula_db, "\\n";
						
						
						my \$query = new Spiders::DatabaseQuery();
						\$query->setBin("$lib");										
						\$query->setDatabase(\$params->{structure_database});
						\$query->setFormula(\$formula_db);
						\$query->setMass(\$mass_db);
						\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
						#\$query->setDatatype("SMILES");
						\$query->setDBname(\$params->{database});						
						\$str_db_return = \$query->QueryStructureDatabase();				
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;	
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Target","\\t",\$formula_db,"\\n";					
						
						if(scalar(\@\$str_db_return)>0)
						{
							foreach (\@\$str_db_return)
							{
								chomp \$_;							
								my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
								print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
							}
						}
						
						\$index++;
					}						
				}
	################################################

				
				
	###### generate decoys 
				my \$formula_unlabel_decoy = "";
				foreach my \$theoretical_formula (\@\$return_decoy) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
						my \$max_C_num = carbon_max_number(\$int_ratio);
											
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
												
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/C\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{

						my \$max_C_num = carbon_max_number(\$int_ratio);						
						\$nitrogen_match = 1;
						for(\$C_i=\$C_num;\$C_i<\$max_C_num;\$C_i++)
						{
							\$carbon_match_temp = (\$theoretical_formula=~\/C\$C_i\\D+\/);
							if(\$carbon_match_temp)
							{
								\$carbon_match = 1;
							}							
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{				
						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
						\$formula_unlabel_decoy .= "\${formula_db},";	

					}	
				}
				
				if(length(\$formula_unlabel_decoy)>1)
				{
					my \$min_mass = \$C_mz - 5*1.007277;
					my \$max_mass = \$C_mz + 5*1.007277;				
					my \$value = \$query->partial_candidate(\$formula_unlabel_decoy, \$min_mass, \$max_mass);
					print \@\$value;
					foreach my \$result (\@\$value)
					{
						my \@data = split(/\\\t/,\$result);
						\$formula_db = \$data[0];
						my \$score = \$data[4];
						\$score=~s/MatchScore\://;			
						next if(\$score<0.9);
					
						my \$updated_formula_db = \$formula_db;
						if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
						{
							my \$Hminus = \$2 - \$decoy_H_num;
							\$updated_formula_db = \$1 . \$Hminus . \$3
						}

									
						my \$updated_mass_db = \$mass_db - 1.00782503207 * \$decoy_H_num;
					
						my \$query = new Spiders::DatabaseQuery();
						\$query->setBin("$lib");
						\$query->setDatabase(\$params->{structure_database});
						\$query->setFormula(\$updated_formula_db);
						\$query->setMass(\$updated_mass_db);
						\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
						#\$query->setDatatype("SMILES");						
						\$query->setDBname(\$params->{database});
						
						my \$str_db_return = \$query->QueryStructureDatabase();				
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;	
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Decoy","\\t",\$formula_db,"\\n";
						
						if(scalar(\@\$str_db_return)>0)
						{
							foreach (\@\$str_db_return)
							{
								chomp \$_;
								my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
								print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
								
								my \$smiles_orig = \$smiles;
								\$index++;

							}
						}						
						
						\$index++;		
					}						
				}																
			}				
		}	
    }
	
}
	
sub carbon_max_number
{
	my (\$int_ratio) = \@_;
	my \$max_C_num = int(\$int_ratio + \$int_ratio*0.1)+1;		
	return \$max_C_num;
}	
	
EOF
	close(RUNSHELL);	
}



sub create_mscript_unlabel
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();

	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;

GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,		
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

if(\$database=~/\,/)
{
	\$decoy_database=~s/\,/\_decoy\,/;
}
my \$falsified_shift = \$params->{'mass_shift'};

my \@dtafiles = \@ARGV;
my \$decoy_H_num = 1;
if(defined(\$params->{'decoy_strategy'}))
{
	\$decoy_H_num = \$params->{'decoy_strategy'} - 1;
}
my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;
    my \$missile_index=1;

    my \$dtafile = abs_path(\$_);
	
	print "\\nProcessing scan: \$dtafile\\n";		
	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\\\/)\(\\d\+\)\\\.MS1\/\$2\/;
#	next if (\$scan ne \$params->{"unlabel_scan"});

	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	
   # print RESULT "MS1 scan: \$dtafile\\n";
  #  print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tFormula\\tPscore\\tStructure (SMILES)\\tType\\n";

	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_diff_score","\\t", "C12_C13_int_rel","\\t", "C12_C13_int_sim","\\t","C12_N15_diff","\\t", "C12_N15_int_rel","\\t", "C12_N15_int_sim", "\\t","Scan number", "\\t","C12 intensity", "\\t","N15 intensity", "\\t","C13 intensity","\\tType\\tFormula\\tAdduct\\n";
    #print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tFormula\\tPscore\\tStructure (SMILES)\\tType\\n";

	print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tPscore\\tFormula\\tAdduct\\tName\\tStructure (SMILES)\\tInChI\\tType\\n";
  
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	my (\$mz_hash_ref,\$int_ratio_hash) = \$deisotope->MS1_deisotope();
	\%mz_hash = \%\$mz_hash_ref;
	
    print scalar(keys \%mz_hash)," features were detected\\n";	
	
    foreach my \$mono_mass (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
#		next if (\$mono_mass < \$params->{"unlabel_min_mass"} or \$mono_mass > \$params->{"unlabel_max_mass"});
		
		my \$int_ratio = \$int_ratio_hash->{\$mono_mass};
        print "\nperforming mass formula database search for mass: \$mono_mass\\n";
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass + 1.007277 + \$shift;
		}
		else
		{
			\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		}
		
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");		
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass - 1.007277;
		}
		else
		{
			\$mono_mass = \$mono_mass + 1.007277;
		}

		\$mono_mass = \$mono_mass + \$falsified_shift;
		
        print "found ",scalar(\@\$return)," formulas, including ",scalar(\@\$return_norm)," targets and ",scalar(\@\$return_decoy)," decoys\\n";	
		
		foreach my \$theoretical_formula (\@\$return) 
		{				
			chomp \$theoretical_formula;

			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
#			my \$carbon_match = carbon_number_matching(\$theoretical_formula,\$int_ratio);
#			next if(\$carbon_match == -1);
			\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
			\$formula_db=~s/(\\D+)1\$/\$1/g;				
			print MRESULT \$index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","pass","\\t",\$formula_db,"\t","NA","\\n";
				
			\$index++;					
		}
					
		foreach my \$theoretical_formula (\@\$return_norm) 
		{				
			chomp \$theoretical_formula;
			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
#			my \$carbon_match = carbon_number_matching(\$theoretical_formula,\$int_ratio);
#			next if(\$carbon_match == -1);
			\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
			\$formula_db=~s/(\\D+)1\$/\$1/g;				
			print MRESULT \$missile_index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","Target","\\t",\$formula_db,"\t","NA","\\n";
			my \$query = new Spiders::DatabaseQuery();
			\$query->setBin("$lib");								

			\$query->setDatabase(\$params->{structure_database});
			\$query->setFormula(\$formula_db);
			\$query->setMass(\$mass_db);
			\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
			\$query->setDBname(\$params->{database});				
			\$str_db_return = \$query->QueryStructureDatabase();		   
			
		   
			if(scalar(\@\$str_db_return)>0)
			{
				foreach (\@\$str_db_return)
				{									
					chomp \$_;
					my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
							
					print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
													
					\$index++;
				}
			}
			\$missile_index++;
			
		}


		foreach my \$theoretical_formula (\@\$return_decoy) 
		{				
			chomp \$theoretical_formula;
			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
#			my \$carbon_match = carbon_number_matching(\$theoretical_formula,\$int_ratio);
#			next if(\$carbon_match == -1);
			\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
			\$formula_db=~s/(\\D+)1\$/\$1/g;				
			print MRESULT \$missile_index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","Decoy","\\t",\$formula_db,"\t","NA","\\n";
			\$index++;

			my \$updated_formula_db = \$formula_db;
			if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
			{
				my \$Hminus = \$2 - \$decoy_H_num;
				\$updated_formula_db = \$1 . \$Hminus . \$3
			}

						
			my \$updated_mass_db = \$mass_db - 1.00782503207 * \$decoy_H_num;
		
			my \$query = new Spiders::DatabaseQuery();
			\$query->setBin("$lib");
			\$query->setDatabase(\$params->{structure_database});
			\$query->setFormula(\$updated_formula_db);
			\$query->setMass(\$updated_mass_db);
			\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
			\$query->setDBname(\$params->{database});				
			\$str_db_return = \$query->QueryStructureDatabase();

			if(scalar(\@\$str_db_return)>0)
			{
				foreach (\@\$str_db_return)
				{
					chomp \$_;
					my (\$formula_db_1,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
					print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tDecoy","\\n";																						
					\$index++;				
				}
			}		
		}
	}				
}

sub carbon_number_matching
{
	my (\$theoretical_formula,\$int_ratio) = \@_;
	my \$carbon_match = -1;			
	if(\$int_ratio != 0)
	{
		my \$min_C_num = int(\$int_ratio - \$int_ratio*0.3);
		if(\$min_C_num < 0)
		{
			\$min_C_num = 0;
		}
		my \$max_C_num = int(\$int_ratio + \$int_ratio*0.3)+1;
		for(\$i=\$min_C_num;\$i<\$max_C_num;\$i++)
		{
			if(\$theoretical_formula=~\/C\$i\\D+\/)
			{
				\$carbon_match = 1;
			}				
		}
	}
	else
	{
		\$carbon_match = 0;
	}
	return \$carbon_match;
}
			
EOF
	close(RUNSHELL);
}	

sub create_mscript_unlabel_library
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();

	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;

GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,		
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

if(\$database=~/\,/)
{
	\$decoy_database=~s/\,/\_decoy\,/;
}
my \$falsified_shift = \$params->{'mass_shift'};

my \@dtafiles = \@ARGV;

my \%libraryhash;
open(LIB,"$lib\/SpectralLib\/library.txt") || die "can not open the library\\n";
while(<LIB>)
{
	my \@data = split(/\\t/,\$_);
	my \$formula = \$data[1];
	\$formula=~s/(\\D+)1(\\D+)/\$1\$2/g;
	\$formula=~s/(\\D+)1\$/\$1/g;	
	\$libraryhash{\$formula}{\$data[5]}=\$_;
}
my \$decoy_H_num = 1;
if(defined(\$params->{'decoy_strategy'}))
{
	\$decoy_H_num = \$params->{'decoy_strategy'} - 1;
}
my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;
    my \$missile_index=1;

    my \$dtafile = abs_path(\$_);
	
	print "\\nProcessing scan: \$dtafile\\n";		
	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\\\/)\(\\d\+\)\\\.MS1\/\$2\/;
#	next if (\$scan ne \$params->{"unlabel_scan"});

	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	

	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_diff_score","\\t", "C12_C13_int_rel","\\t", "C12_C13_int_sim","\\t","C12_N15_diff","\\t", "C12_N15_int_rel","\\t", "C12_N15_int_sim", "\\t","Scan number", "\\t","C12 intensity", "\\t","N15 intensity", "\\t","C13 intensity","\\tType\\tFormula\\tAdduct\\n";

	print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tPscore\\tFormula\\tAdduct\\tName\\tStructure (SMILES)\\tInChI\\tType\\n";
  
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	my (\$mz_hash_ref,\$int_ratio_hash) = \$deisotope->MS1_deisotope();
	\%mz_hash = \%\$mz_hash_ref;
	
    print scalar(keys \%mz_hash)," features were detected\\n";	
	
    foreach my \$mono_mass (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
#		next if (\$mono_mass < \$params->{"unlabel_min_mass"} or \$mono_mass > \$params->{"unlabel_max_mass"});
		
		my \$int_ratio = \$int_ratio_hash->{\$mono_mass};
        print "\nperforming mass formula database search for mass: \$mono_mass\\n";
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass + 1.007277 + \$shift;
		}
		else
		{
			\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		}
		
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");		
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass - 1.007277;
		}
		else
		{
			\$mono_mass = \$mono_mass + 1.007277;
		}

		\$mono_mass = \$mono_mass + \$falsified_shift;
		
        print "found ",scalar(\@\$return)," formulas, including ",scalar(\@\$return_norm)," targets and ",scalar(\@\$return_decoy)," decoys\\n";	
		
		foreach my \$theoretical_formula (\@\$return) 
		{				
			chomp \$theoretical_formula;

			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
			my \$carbon_match = carbon_number_matching(\$theoretical_formula,\$int_ratio);
			
			next if(\$carbon_match == -1);
			\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
			\$formula_db=~s/(\\D+)1\$/\$1/g;			
			next if(!defined(\$libraryhash{\$formula_db}));			
			print MRESULT \$index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","pass","\\t",\$formula_db,"\\t","NA","\\n";
				
			\$index++;					
		}
					
		foreach my \$theoretical_formula (\@\$return_norm) 
		{				
			chomp \$theoretical_formula;
			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
			my \$carbon_match = carbon_number_matching(\$theoretical_formula,\$int_ratio);
			next if(\$carbon_match == -1);
			\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
			\$formula_db=~s/(\\D+)1\$/\$1/g;			
			next if(!defined(\$libraryhash{\$formula_db}));			
			print MRESULT \$missile_index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","Target","\\t",\$formula_db,"\\t","NA","\\n";
			my \$query = new Spiders::DatabaseQuery();
			\$query->setBin("$lib");								

			\$query->setDatabase(\$params->{structure_database});
			\$query->setFormula(\$formula_db);
			\$query->setMass(\$mass_db);
			\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
			\$query->setDBname(\$params->{database});				
			#\$str_db_return = \$query->QueryStructureDatabase();		   
			
		   
			foreach my \$smiles (keys \%{\$libraryhash{\$formula_db}})
			{
				my \$line = \$libraryhash{\$formula_db}{\$smiles};
				chop \$line;
				my \@data = split(/\\t/,\$line);
				my (\$InChI,\$IUPAC)=(\$data[3],\$data[6]);				
				print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";																						
				\$index++;				
			}	
			\$missile_index++;
			
		}


		foreach my \$theoretical_formula (\@\$return_decoy) 
		{				
			chomp \$theoretical_formula;
			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
			my \$carbon_match = carbon_number_matching(\$theoretical_formula,\$int_ratio);
			next if(\$carbon_match == -1);
			\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
			\$formula_db=~s/(\\D+)1\$/\$1/g;			
			next if(!defined(\$libraryhash{\$formula_db}));			
			print MRESULT \$missile_index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","Decoy","\\t",\$formula_db,"\\t","NA","\\n";
			\$index++;

			my \$updated_formula_db = \$formula_db;
			if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
			{
				my \$Hminus = \$2 - \$decoy_H_num;
				\$updated_formula_db = \$1 . \$Hminus . \$3
			}

						
			my \$updated_mass_db = \$mass_db - 1.00782503207 * \$decoy_H_num;
		
			my \$query = new Spiders::DatabaseQuery();
			\$query->setBin("$lib");
			\$query->setDatabase(\$params->{structure_database});
			\$query->setFormula(\$updated_formula_db);
			\$query->setMass(\$updated_mass_db);
			\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
			\$query->setDBname(\$params->{database});				
			#\$str_db_return = \$query->QueryStructureDatabase();


			foreach my \$smiles (keys \%{\$libraryhash{\$formula_db}})
			{
				my \$line = \$libraryhash{\$formula_db}{\$smiles};
				chop \$line;
				my \@data = split(/\\t/,\$line);
				my (\$InChI,\$IUPAC)=(\$data[3],\$data[6]);				
				print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";																						
				\$index++;				
			}		
		}
	}				
}

sub carbon_number_matching
{
	my (\$theoretical_formula,\$int_ratio) = \@_;
	my \$carbon_match = -1;			
	if(\$int_ratio != 0)
	{
		my \$min_C_num = int(\$int_ratio - \$int_ratio*0.3);
		if(\$min_C_num < 0)
		{
			\$min_C_num = 0;
		}
		my \$max_C_num = int(\$int_ratio + \$int_ratio*0.3)+1;
		for(\$i=\$min_C_num;\$i<\$max_C_num;\$i++)
		{
			if(\$theoretical_formula=~\/C\$i\\D+\/)
			{
				\$carbon_match = 1;
			}				
		}
	}
	else
	{
		\$carbon_match = 0;
	}
	return \$carbon_match;
}
			
EOF
	close(RUNSHELL);
}	

sub create_mscript_unlabel_dta
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();

	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;

GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,		
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

if(\$database=~/\,/)
{
	\$decoy_database=~s/\,/\_decoy\,/;
}
my \$falsified_shift = \$params->{'mass_shift'};

my \@dtafiles = \@ARGV;

my \$ms_hash_mol;
my \$decoy_H_num = 1;
if(defined(\$params->{'decoy_strategy'}))
{
	\$decoy_H_num = \$params->{'decoy_strategy'} - 1;
}
foreach(\@dtafiles)
{
    my \$index=1;
    my \$missile_index=1;

    my \$dtafile = abs_path(\$_);
	
	print "\\nProcessing scan: \$dtafile\\n";		
	my \$dta=new Spiders::Dta;  
	\$dta->set_dta_file(\$dtafile) ;
	\$dta->process_dtafile(\$dtafile);
	\$mono_mass = \$dta->get_prec_mz()-1.007277;
    my \$scan = basename(\$dtafile);
    	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	
	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_diff_score","\\t", "C12_C13_int_rel","\\t", "C12_C13_int_sim","\\t","C12_N15_diff","\\t", "C12_N15_int_rel","\\t", "C12_N15_int_sim", "\\t","Scan number", "\\t","C12 intensity", "\\t","N15 intensity", "\\t","C13 intensity","\\tType\\tFormula\\tAdduct\\n";

	print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tPscore\\tFormula\\tAdduct\\tName\\tStructure (SMILES)\\tInChI\\tType\\n";
  
	
    print "\\nperforming mass formula database search for mass: \$mono_mass\\n";
		
	my \$query = new Spiders::DatabaseQuery();
	\$query->setBin("$lib");		
	my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);
	
	
    print "found ",scalar(\@\$return)," formulas, including ",scalar(\@\$return_norm)," targets and ",scalar(\@\$return_decoy)," decoys\\n";	
	
	foreach my \$theoretical_formula (\@\$return) 
	{				
		chomp \$theoretical_formula;

		my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

		\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
		\$formula_db=~s/(\\D+)1\$/\$1/g;					
		print MRESULT \$index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","pass","\\t",\$formula_db,"\\t","NA","\\n";
			
		\$index++;					
	}
				

					
	foreach my \$theoretical_formula (\@\$return_norm) 
	{				
		chomp \$theoretical_formula;
		my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

		\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
		\$formula_db=~s/(\\D+)1\$/\$1/g;				
		print MRESULT \$missile_index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","Target","\\t",\$formula_db,"\t","NA","\\n";
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");								

		\$query->setDatabase(\$params->{structure_database});
		\$query->setFormula(\$formula_db);
		\$query->setMass(\$mass_db);
		\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
		\$query->setDBname(\$params->{database});				
		\$str_db_return = \$query->QueryStructureDatabase();		   
		
	   
		if(scalar(\@\$str_db_return)>0)
		{
			foreach (\@\$str_db_return)
			{									
				chomp \$_;
				my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
						
				print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
												
				\$index++;
			}
		}
		\$missile_index++;
		
	}


	foreach my \$theoretical_formula (\@\$return_decoy) 
	{				
		chomp \$theoretical_formula;
		my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

		\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
		\$formula_db=~s/(\\D+)1\$/\$1/g;				
		print MRESULT \$missile_index,"\\t",\$mono_mass,"\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t","NA","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","NA","\\t","NA","\\t","Decoy","\\t",\$formula_db,"\t","NA","\\n";
		\$index++;

		my \$updated_formula_db = \$formula_db;
		if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
		{
			my \$Hminus = \$2 - 1* \$decoy_H_num;
			\$updated_formula_db = \$1 . \$Hminus . \$3
		}

					
		my \$updated_mass_db = \$mass_db - 1.00782503207 * \$decoy_H_num;
	
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");
		\$query->setDatabase(\$params->{structure_database});
		\$query->setFormula(\$updated_formula_db);
		\$query->setMass(\$updated_mass_db);
		\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
		\$query->setDBname(\$params->{database});				
		\$str_db_return = \$query->QueryStructureDatabase();

		if(scalar(\@\$str_db_return)>0)
		{
			foreach (\@\$str_db_return)
			{
				chomp \$_;
				my (\$formula_db_1,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
				print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$mz_hash{\$mono_mass},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tDecoy","\\n";																						
				\$index++;				
			}
		}		
	}
}
		
EOF
	close(RUNSHELL);
}	

sub create_mscript_library
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	my $R_lib = $self->get_R_library_path();
	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,
		'-NC=f'=>\\\$NC,
		'-CC=f'=>\\\$CC,
		'-NC_std=f'=>\\\$NC_std,
		'-CC_std=f'=>\\\$CC_std,
		'-NC_defect_loc=f'=>\\\$NC_defect_loc,
		'-CC_defect_loc=f'=>\\\$CC_defect_loc,
		'-NC_defect_scale=f'=>\\\$NC_defect_scale,
		'-CC_defect_scale=f'=>\\\$CC_defect_scale,
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

if(\$database=~/\,/)
{
	\$decoy_database=~s/\,/\_decoy\,/;
}
if(\$params->{'mass_correction}'}==1)
{
	my \$mass_shift = retrieve("$dir/.mass_shift");
}

my \$falsified_shift = \$params->{'mass_shift'};

my \@dtafiles = \@ARGV;
my \$decoy_H_num = 1;
if(defined(\$params->{'decoy_strategy'}))
{
	\$decoy_H_num = \$params->{'decoy_strategy'} - 1;
}
my \$ms_hash_mol;

my \%libraryhash;
open(LIB,"\$lib\/SpectralLib\/library.txt") || die "can not open the library\\n";
while(<LIB>)
{
	my \@data = split(/\\t/,\$_);
	my \$formula = \$data[1];
	\$formula=~s/(\\D+)1(\\D+)/\$1\$2/g;
	\$formula=~s/(\\D+)1\$/\$1/g;	
	\$libraryhash{\$formula}{\$data[5]}=\$_;
}

foreach(\@dtafiles)
{
    my \$index=1;
    my \$dtafile = abs_path(\$_);
	print "\\nProcessing scan: \$dtafile\\n";	
	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\)\\\.\(\\d\+\)\\\.dta\/\$2\/;	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_diff_score","\\t", "C12_C13_int_rel","\\t", "C12_C13_int_sim","\\t","C12_N15_diff","\\t", "C12_N15_int_rel","\\t", "C12_N15_int_sim", "\\t","Scan number", "\\t","C12 intensity","\\t","N15 intensity","\\t","C13 intensity","\\tType\\tFormula\\tAdduct\\n";
		
   # print RESULT "MS1 scan: \$dtafile\\n";
    print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tPscore\\tFormula\\tAdduct\\tName\\tStructure (SMILES)\\tInChI\\tType\\n";
	
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	my (\$mz_hash_ref,\$int_ratio_hash,\$peak_C,\$charge_hash) = \$deisotope->MS1_deisotope();
    #my \$dta=new Spiders::Dta;  
    #\$dta->set_dta_file(\$dtafile) ;
    #\$dta->process_dtafile(\$dtafile);

    #my \%mz_hash = %{\$dta->{'_mz_int_hash'}};
	\%mz_hash = \%\$mz_hash_ref;	
    print scalar(keys \%mz_hash)," features were detected\\n";		
	
	print "Peaks pairing\\n";			
    my \$cluster;
	my \$pair = new Spiders::Pairing();
	\$pair->set_parameter(\$params);
	\$pair->set_library_path("$R_lib");
	\$pair->set_NC_ratio(\$NC);
	\$pair->set_CC_ratio(\$CC);
	\$pair->set_NC_std(\$NC_std);
	\$pair->set_CC_std(\$CC_std);
	\$pair->set_NC_defect_loc(\$NC_defect_loc);
	\$pair->set_CC_defect_loc(\$CC_defect_loc);
	\$pair->set_NC_defect_scale(\$NC_defect_scale);
	\$pair->set_CC_defect_scale(\$CC_defect_scale);
	
    foreach my \$mz (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
        next if (\!defined(\$mz_hash{\$mz}));
		\$cluster->{\$mz}->{0}->{'mz'} = \$mz;
        \$cluster->{\$mz}->{0}->{'int'} = \$mz_hash{\$mz};

        \$pair->clustering(\$cluster,\\\%mz_hash,\$mz);
    }
	
	
############### only for summary purpose ########
    my \$dta=new Spiders::Dta;                  #
	#\$dtafile .= ".iso";                        #
	\$dta->set_prec_mz("1");                    #
    \$dta->set_charge("0");	                    #
    \$dta->set_dta_file(\$dtafile) ;	        #
	\$dta->set_mz_int_hash(\\\%mz_hash);	
	\$dta->write_dta_file();                    #
#################################################
	
	
	
	my \$dirname  = dirname(\$dtafile);
	
    my \$cluster_mono = \$pair->select_mono_peaks(\$cluster);

    my (\$formula_comb) = \$pair->pairing(\$cluster_mono,\$dirname,\$dtafile,\$peak_C,\$charge_hash);
	
	\$index=1;	
	
	print scalar(keys \%\$formula_comb)," MISSILES were detected\\n";
    foreach my \$mono_mass (sort {\$a<=>\$b} keys \%\$formula_comb)
    {
	
		my \$int_ratio = \$int_ratio_hash->{\$mono_mass};
        print "\nperforming mass formula database search for mass: \$mono_mass\\n";
		
        my \@missile_formula = \@{\$formula_comb->{\$mono_mass}->{'number'}};
        next if(\$mono_mass>1500);
		next if(\$#missile_formula<0);
###### For decoy testing		
		my \$shift = 0;
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass + 1.007277 + \$shift;
		}
		else
		{
			\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		}
		if(\$params->{mass_correction})
		{
			\$mono_mass = \$mono_mass - \$mono_mass * \$mass_shift->{\$scan} / 1000000;
		}		
		\$mono_mass = \$mono_mass + \$falsified_shift;
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");		
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);		
		# my \$return =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{mass_formula_database}, "yes");
		# my \$return_decoy =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{decoy_str_mass_formula_database}, "all");
		# my \$return_norm =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all");
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass - 1.007277;
		}
		else
		{
			\$mono_mass = \$mono_mass + 1.007277;
		}
		
##### if the case of "mass only" 
		if(\$params->{'labeled_ID_method'} == 1)
		{
			foreach my \$theoretical_formula (\@\$return) 
			{				
				chomp \$theoretical_formula;
				
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				next if(!defined(\$libraryhash{\$formula_db}));
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\n";						
				\$index++;					
			}
									
			foreach my \$theoretical_formula (\@\$return_norm) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				next if(!defined(\$libraryhash{\$formula_db}));				
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\tnorm\\n";

									
				\$query->setDatabase(\$params->{structure_database});
				\$query->setFormula(\$formula_db);
				\$query->setMass(\$mass_db);
				#\$query->setDatatype("ID,INCHIKEY,SMILES,GENERALNAME");

				\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
				#\$query->setDatatype("SMILES");
				\$query->setDBname(\$params->{database});				
				\$str_db_return = \$query->QueryStructureDatabase();				
	
				if(scalar(\@\$str_db_return)>0)
				{
					foreach (\@\$str_db_return)
					{
						chomp \$_;					
						my (\$ID,\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
						next if(!defined(\$libraryhash{\$formula_db}{\$smiles}));
						print RESULT \$ID,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$formula_db,"\\t",  "","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
						
#						print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$formula_db,"\\t",  "","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
					}
				}
						
				\$index++;					
			}


			foreach my \$theoretical_formula (\@\$return_decoy) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				next if(!defined(\$libraryhash{\$formula_db}));				
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\tdecoy\\n";		
				\$index++;					
			}			
			
			
		}
		else
		{
			my \$detected_formula=0;
			my \$adduct_formula=0;			
			foreach my \$missile_formula (\@missile_formula)
			{
				my \@N_C = split(/\\,/,\$missile_formula);
				my (\$C_num_mz,\$c_diff_score,\$c_int_score, \$c_similarity) = split(/\\|/,\$N_C[0]);
				my (\$N_num_mz,\$n_diff_score,\$n_int_score, \$n_similarity) = split(/\\|/,\$N_C[1]);
				my \$pair_score = \$N_C[2];
				my (\$C_num,\$C_mz) = split(/\\:/,\$C_num_mz);
				my (\$N_num,\$N_mz) = split(/\\:/,\$N_num_mz);	

###### correct mass before the search 		
				if(\$params->{mass_correction})
				{
					\$C_mz = \$C_mz - \$C_mz * \$mass_shift->{\$scan} / 1000000;
					\$N_mz = \$N_mz - \$N_mz * \$mass_shift->{\$scan} / 1000000;					
				}
				
					
				foreach my \$theoretical_formula (\@\$return) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = 1;
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\\n";
						exit(1);
					}					
										
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{
						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;						
						next if(!defined(\$libraryhash{\$formula_db}));
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t", "pass", "\\t", \$formula_db, "\\t \\n";
						\$detected_formula++;						
					}
								
				}
				if(\$detected_formula==0 and \$params->{'adduct'} == 1)
				{
					my \$adduct_hash = \$p->get_adducts(\$params);
					foreach \$adduct_name (keys \%\$adduct_hash)
					{
						my \$mono_mass_no_adduct = \$mono_mass - \$adduct_hash->{\$adduct_name};
						my \$N_mz_no_adduct = \$C_mz_no_adduct = 0;
						if(\$N_mz>0)
						{
							\$N_mz_no_adduct = \$N_mz - \$adduct_hash->{\$adduct_name};
						}
						if(\$C_mz>0)
						{
							\$C_mz_no_adduct = \$C_mz - \$adduct_hash->{\$adduct_name};
						}					
						if(\$params->{'mode'} == -1)
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct + 1.007277;
						}
						else
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct - 1.007277;
						}					
						(\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass_no_adduct, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);	
						if(\$params->{'mode'} == -1)
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct - 1.007277;
						}
						else
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct + 1.007277;
						}						
						
						foreach my \$theoretical_formula (\@\$return) 
						{
							chomp \$theoretical_formula;
							my \$carbon_match = 0;
							my \$nitrogen_match = 0;
							
							if(\$params->{'labeled_ID_method'} == 2)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
								if(\$N_num eq "N0")
								{
									if(\$theoretical_formula!~/N/)
									{
										\$nitrogen_match=1;
									}						
								}
							}
							elsif(\$params->{'labeled_ID_method'} == 3)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = 1;
							}
							elsif(\$params->{'labeled_ID_method'} == 4)
							{
								\$carbon_match = 1;
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\\\n";
								exit(1);
							}					
												
							if(\$carbon_match == 1 and \$nitrogen_match ==1)
							{
								
								my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
								my \$adduct_mass = \$adduct_hash->{\$adduct_name};
								next if(!defined(\$libraryhash{\$formula_db}));
								print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t", "pass", "\\t", \$formula_db, "\\t",  \$adduct_name,"\\n";
								\$adduct_formula++;						
							}
										
						}
						
						foreach my \$theoretical_formula (\@\$return_norm) 
						{
						
						
							chomp \$theoretical_formula;				
							my \$carbon_match = 0;
							my \$nitrogen_match = 0;
							
							if(\$params->{'labeled_ID_method'} == 2)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
								if(\$N_num eq "N0")
								{
									if(\$theoretical_formula!~/N/)
									{
										\$nitrogen_match=1;
									}						
								}
							}
							elsif(\$params->{'labeled_ID_method'} == 3)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = 1;
							}
							elsif(\$params->{'labeled_ID_method'} == 4)
							{
								\$carbon_match = 1;
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\\n";
								exit(1);
							}					
								
							if(\$carbon_match == 1 and \$nitrogen_match ==1)
							{
							
								my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

								my \$query = new Spiders::DatabaseQuery();
								\$query->setBin("$lib");										
								\$query->setDatabase(\$params->{structure_database});
								\$query->setFormula(\$formula_db);
								\$query->setMass(\$mass_db);
								\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
								#\$query->setDatatype("SMILES");
								\$query->setDBname(\$params->{database});						
								\$str_db_return = \$query->QueryStructureDatabase();				
								my \$adduct = (\$adduct_formula > 0) ? \$adduct_name : "";
								my \$adduct_mass = (\$adduct_formula > 0) ? \$adduct_hash->{\$adduct_name} : 0;
								\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
								\$formula_db=~s/(\\D+)1\$/\$1/g;
								next if(!defined(\$libraryhash{\$formula_db}));								
								print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Target","\\t",\$formula_db, "\\t",  \$adduct,"\\n";					
								
								if(scalar(\@\$str_db_return)>0)
								{
									foreach (\@\$str_db_return)
									{
										chomp \$_;							
										my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);

										next if(!defined(\$libraryhash{\$formula_db}{\$smiles}));										
										print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  \$adduct,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
										
										
									}
								}
								
								\$index++;		
							}	
						}
			################################################

						
						
			###### generate decoys 
						foreach my \$theoretical_formula (\@\$return_decoy) 
						{
							chomp \$theoretical_formula;
							my \$carbon_match = 0;
							my \$nitrogen_match = 0;
							
							if(\$params->{'labeled_ID_method'} == 2)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
								if(\$N_num eq "N0")
								{
									if(\$theoretical_formula!~/N/)
									{
										\$nitrogen_match=1;
									}						
								}
							}
							elsif(\$params->{'labeled_ID_method'} == 3)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = 1;
							}
							elsif(\$params->{'labeled_ID_method'} == 4)
							{
								\$carbon_match = 1;
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\\n";
								exit(1);
							}					
								
							if(\$carbon_match == 1 and \$nitrogen_match ==1)
							{				
								
								my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

								my \$updated_formula_db = \$formula_db;
								if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
								{
									my \$Hminus = \$2 - 1 * \$decoy_H_num;
									\$updated_formula_db = \$1 . \$Hminus . \$3
								}

											
								my \$updated_mass_db = \$mass_db - 1.00782503207  * \$decoy_H_num;
									
								my \$query = new Spiders::DatabaseQuery();
								\$query->setBin("$lib");
								\$query->setDatabase(\$params->{structure_database});
								\$query->setFormula(\$updated_formula_db);
								\$query->setMass(\$updated_mass_db);
								\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
								#\$query->setDatatype("SMILES");						
								\$query->setDBname(\$params->{database});
								
								my \$str_db_return = \$query->QueryStructureDatabase();	
								my \$adduct = (\$adduct_formula > 0) ? \$adduct_name : "";
								\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
								\$formula_db=~s/(\\D+)1\$/\$1/g;									
								my \$adduct_mass = (\$adduct_formula > 0) ? \$adduct_hash->{\$adduct_name} : 0; 

								next if(!defined(\$libraryhash{\$formula_db}));									
								print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Decoy","\\t",\$formula_db, "\\t",  \$adduct,"\\n";
								
								if(scalar(\@\$str_db_return)>0)
								{
									foreach (\@\$str_db_return)
									{
										chomp \$_;
										my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
										next if(!defined(\$libraryhash{\$formula_db}{\$smiles}));
										print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  \$adduct,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tDecoy","\\n";

										my \$smiles_orig = \$smiles;
										\$index++;

									}
								}						
								
								\$index++;		
							}	
						}
					}
	
				}
				else
				{

####### generating structure and scoring #######################
			
					foreach my \$theoretical_formula (\@\$return_norm) 
					{
					
					
						chomp \$theoretical_formula;				
						my \$carbon_match = 0;
						my \$nitrogen_match = 0;
						
						if(\$params->{'labeled_ID_method'} == 2)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							if(\$N_num eq "N0")
							{
								if(\$theoretical_formula!~/N/)
								{
									\$nitrogen_match=1;
								}						
							}
						}
						elsif(\$params->{'labeled_ID_method'} == 3)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = 1;
						}
						elsif(\$params->{'labeled_ID_method'} == 4)
						{
							\$carbon_match = 1;
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
						}
						else
						{
							print "please set the labeled_ID_method parameter [1-4]\\n";
							exit(1);
						}					
							
						if(\$carbon_match == 1 and \$nitrogen_match ==1)
						{
						
							my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

							my \$query = new Spiders::DatabaseQuery();
							\$query->setBin("$lib");										
							\$query->setDatabase(\$params->{structure_database});
							\$query->setFormula(\$formula_db);
							\$query->setMass(\$mass_db);
							\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
							#\$query->setDatatype("SMILES");
							\$query->setDBname(\$params->{database});						
							\$str_db_return = \$query->QueryStructureDatabase();				

							\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
							\$formula_db=~s/(\\D+)1\$/\$1/g;

							next if(!defined(\$libraryhash{\$formula_db}));								
							print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Target","\\t",\$formula_db, "\\t",  \$adduct,"\\n";					
							
							if(scalar(\@\$str_db_return)>0)
							{
								foreach (\@\$str_db_return)
								{
									chomp \$_;							
									my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
									next if(!defined(\$libraryhash{\$formula_db}{\$smiles}));
									print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  "","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
									
								}
							}
							
							\$index++;		
						}	
					}
		################################################

					
					
		###### generate decoys 
					foreach my \$theoretical_formula (\@\$return_decoy) 
					{
						chomp \$theoretical_formula;
						my \$carbon_match = 0;
						my \$nitrogen_match = 0;
						
						if(\$params->{'labeled_ID_method'} == 2)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							if(\$N_num eq "N0")
							{
								if(\$theoretical_formula!~/N/)
								{
									\$nitrogen_match=1;
								}						
							}
						}
						elsif(\$params->{'labeled_ID_method'} == 3)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = 1;
						}
						elsif(\$params->{'labeled_ID_method'} == 4)
						{
							\$carbon_match = 1;
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
						}
						else
						{
							print "please set the labeled_ID_method parameter [1-4]\\n";
							exit(1);
						}					
							
						if(\$carbon_match == 1 and \$nitrogen_match ==1)
						{				
							
							my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

							my \$updated_formula_db = \$formula_db;
							if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
							{
								my \$Hminus = \$2 - 1  * \$decoy_H_num;
								\$updated_formula_db = \$1 . \$Hminus . \$3
							}

										
							my \$updated_mass_db = \$mass_db - 1.00782503207 * \$decoy_H_num;

							my \$query = new Spiders::DatabaseQuery();
							\$query->setBin("$lib");
							\$query->setDatabase(\$params->{structure_database});
							\$query->setFormula(\$updated_formula_db);
							\$query->setMass(\$updated_mass_db);
							\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
							#\$query->setDatatype("SMILES");						
							\$query->setDBname(\$params->{database});
							
							my \$str_db_return = \$query->QueryStructureDatabase();	
							\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
							\$formula_db=~s/(\\D+)1\$/\$1/g;	

							next if(!defined(\$libraryhash{\$formula_db}));								
							print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Decoy","\\t",\$formula_db, "\\t",  \$adduct,"\\n";
							
							if(scalar(\@\$str_db_return)>0)
							{
								foreach (\@\$str_db_return)
								{
									chomp \$_;
									my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
									next if(!defined(\$libraryhash{\$formula_db}{\$smiles}));
									print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  "","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tDecoy","\\n";
									my \$smiles_orig = \$smiles;
									\$index++;

								}
							}						
							
							\$index++;		
						}	
					}
				}
			}			
		}	
    }
	
}
		
EOF
	close(RUNSHELL);
}
	
sub create_mscript
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	my $R_lib = $self->get_R_library_path();

	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,
		'-NC=f'=>\\\$NC,
		'-CC=f'=>\\\$CC,
		'-NC_std=f'=>\\\$NC_std,
		'-CC_std=f'=>\\\$CC_std,
		'-NC_defect_loc=f'=>\\\$NC_defect_loc,
		'-CC_defect_loc=f'=>\\\$CC_defect_loc,
		'-NC_defect_scale=f'=>\\\$NC_defect_scale,
		'-CC_defect_scale=f'=>\\\$CC_defect_scale,
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

if(\$database=~/\,/)
{
	\$decoy_database=~s/\,/\_decoy\,/;
}
if(\$params->{'mass_correction}'}==1)
{
	my \$mass_shift = retrieve("$dir/.mass_shift");
}

my \$falsified_shift = \$params->{'mass_shift'};

my \@dtafiles = \@ARGV;
my \$decoy_H_num = 1;
if(defined(\$params->{'decoy_strategy'}))
{
	\$decoy_H_num = \$params->{'decoy_strategy'} - 1;
}
my \$ms_hash_mol;

foreach(\@dtafiles)
{
    my \$index=1;
    my \$dtafile = abs_path(\$_);
	print "\\nProcessing scan: \$dtafile\\n";	
	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\)\\\.\(\\d\+\)\\\.dta\/\$2\/;	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_diff_score","\\t", "C12_C13_int_rel","\\t", "C12_C13_int_sim","\\t","C12_N15_diff","\\t", "C12_N15_int_rel","\\t", "C12_N15_int_sim", "\\t","Scan number", "\\t","C12 intensity","\\t","N15 intensity","\\t","C13 intensity","\\tType\\tFormula\\tAdduct\\n";
		
   # print RESULT "MS1 scan: \$dtafile\\n";
    print RESULT "Index\\tC12 mass\\tN15 mass\\tC13 mass\\tMS1 Scan\\tC12 Intensity\\tN15 Intensity\\tC13 Intensity\\tPscore\\tFormula\\tAdduct\\tName\\tStructure (SMILES)\\tInChI\\tType\\n";
	
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	my (\$mz_hash_ref,\$int_ratio_hash,\$peak_C,\$charge_hash) = \$deisotope->MS1_deisotope();
    #my \$dta=new Spiders::Dta;  
    #\$dta->set_dta_file(\$dtafile) ;
    #\$dta->process_dtafile(\$dtafile);

    #my \%mz_hash = %{\$dta->{'_mz_int_hash'}};
	\%mz_hash = \%\$mz_hash_ref;	
    print scalar(keys \%mz_hash)," features were detected\\n";		
	
	print "Peaks pairing\\n";			
    my \$cluster;
	my \$pair = new Spiders::Pairing();
	\$pair->set_parameter(\$params);
	\$pair->set_library_path("$R_lib");
	\$pair->set_NC_ratio(\$NC);
	\$pair->set_CC_ratio(\$CC);
	\$pair->set_NC_std(\$NC_std);
	\$pair->set_CC_std(\$CC_std);
	\$pair->set_NC_defect_loc(\$NC_defect_loc);
	\$pair->set_CC_defect_loc(\$CC_defect_loc);
	\$pair->set_NC_defect_scale(\$NC_defect_scale);
	\$pair->set_CC_defect_scale(\$CC_defect_scale);
	
    foreach my \$mz (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
        next if (\!defined(\$mz_hash{\$mz}));
		\$cluster->{\$mz}->{0}->{'mz'} = \$mz;
        \$cluster->{\$mz}->{0}->{'int'} = \$mz_hash{\$mz};

        \$pair->clustering(\$cluster,\\\%mz_hash,\$mz);
    }
	
	
############### only for summary purpose ########
    my \$dta=new Spiders::Dta;                  #
	#\$dtafile .= ".iso";                        #
	\$dta->set_prec_mz("1");                    #
    \$dta->set_charge("0");	                    #
    \$dta->set_dta_file(\$dtafile) ;	        #
	\$dta->set_mz_int_hash(\\\%mz_hash);	
	\$dta->write_dta_file();                    #
#################################################
	
	
	
	my \$dirname  = dirname(\$dtafile);
	
    my \$cluster_mono = \$pair->select_mono_peaks(\$cluster);

    my (\$formula_comb) = \$pair->pairing(\$cluster_mono,\$dirname,\$dtafile,\$peak_C,\$charge_hash);
	
	\$index=1;	
	
	print scalar(keys \%\$formula_comb)," MISSILES were detected\\n";
    foreach my \$mono_mass (sort {\$a<=>\$b} keys \%\$formula_comb)
    {
	
		my \$int_ratio = \$int_ratio_hash->{\$mono_mass};
        print "\nperforming mass formula database search for mass: \$mono_mass\\n";
		
        my \@missile_formula = \@{\$formula_comb->{\$mono_mass}->{'number'}};
        next if(\$mono_mass>1500);
		next if(\$#missile_formula<0);
###### For decoy testing		
		my \$shift = 0;
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass + 1.007277 + \$shift;
		}
		else
		{
			\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		}
		if(\$params->{mass_correction})
		{
			\$mono_mass = \$mono_mass - \$mono_mass * \$mass_shift->{\$scan} / 1000000;
		}		
		\$mono_mass = \$mono_mass + \$falsified_shift;
		my \$query = new Spiders::DatabaseQuery();
		\$query->setBin("$lib");		
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);		
		# my \$return =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{mass_formula_database}, "yes");
		# my \$return_decoy =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{decoy_str_mass_formula_database}, "all");
		# my \$return_norm =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all");
		
		if(\$params->{'mode'} == -1)
		{
			\$mono_mass = \$mono_mass - 1.007277;
		}
		else
		{
			\$mono_mass = \$mono_mass + 1.007277;
		}
		
##### if the case of "mass only" 
		if(\$params->{'labeled_ID_method'} == 1)
		{
			foreach my \$theoretical_formula (\@\$return) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

				print MRESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t","0","\\t","0","\\t","0","\\t", "0","\\t", "0","\\t","0","\\t", "0","\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t","0","\\t","0","\\t", "pass", "\\t", \$formula_db, "\\t \\n";					
				\$index++;					
			}
									
			foreach my \$theoretical_formula (\@\$return_norm) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				print MRESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t","0","\\t","0","\\t","0","\\t", "0","\\t", "0","\\t","0","\\t", "0","\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t","0","\\t","0","\\t", "Target", "\\t", \$formula_db, "\\t \\n";

									
				\$query->setDatabase(\$params->{structure_database});
				\$query->setFormula(\$formula_db);
				\$query->setMass(\$mass_db);
				#\$query->setDatatype("ID,INCHIKEY,SMILES,GENERALNAME");

				\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
				#\$query->setDatatype("SMILES");
				\$query->setDBname(\$params->{database});				
				\$str_db_return = \$query->QueryStructureDatabase();				
	
				if(scalar(\@\$str_db_return)>0)
				{
					foreach (\@\$str_db_return)
					{
						chomp \$_;					
						my (\$ID,\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
						print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t","0","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";						
					}
				}
						
				\$index++;					
			}


			foreach my \$theoretical_formula (\@\$return_decoy) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
	

				if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
				{
					my \$Hminus = \$2 - 1  * \$decoy_H_num;
					\$updated_formula_db = \$1 . \$Hminus . \$3
				}
	
							
				my \$updated_mass_db = \$mass_db - 1.00782503207  * \$decoy_H_num;				
								
				my \$query = new Spiders::DatabaseQuery();
				\$query->setBin("$lib");
				\$query->setDatabase(\$params->{structure_database});
				\$query->setFormula(\$updated_formula_db);
				\$query->setMass(\$updated_mass_db);
				\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
				#\$query->setDatatype("SMILES");						
				\$query->setDBname(\$params->{database});
				
				my \$str_db_return = \$query->QueryStructureDatabase();	
				my \$adduct = (\$adduct_formula > 0) ? \$adduct_name : "";
				\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
				\$formula_db=~s/(\\D+)1\$/\$1/g;									
				my \$adduct_mass = (\$adduct_formula > 0) ? \$adduct_hash->{\$adduct_name} : 0;                    						
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",0,"\\t",0,"\\t",0,"\\t",0,"\\t",0,"\\t", 0,"\\t", 0,"\\t",0,"\\t", 0,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",0,"\\t",0,"\\t", "Decoy", "\\t", \$formula_db, "\\t \\n";	
								
				if(scalar(\@\$str_db_return)>0)
				{
					foreach (\@\$str_db_return)
					{
						chomp \$_;
						my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);

						print RESULT \$index,"\\t",\$mono_mass,"\\t","0","\\t","0","\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t","0","\\t","0","\\t","0","\\t",\$formula_db,"\\t", "0","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tDecoy","\\n";

						\$index++;

					}
				}									
			}			
						
		}
		else
		{
			my \$detected_formula=0;
			my \$adduct_formula=0;			
			foreach my \$missile_formula (\@missile_formula)
			{
				my \@N_C = split(/\\,/,\$missile_formula);
				my (\$C_num_mz,\$c_diff_score,\$c_int_score, \$c_similarity) = split(/\\|/,\$N_C[0]);
				my (\$N_num_mz,\$n_diff_score,\$n_int_score, \$n_similarity) = split(/\\|/,\$N_C[1]);
				my \$pair_score = \$N_C[2];
				my (\$C_num,\$C_mz) = split(/\\:/,\$C_num_mz);
				my (\$N_num,\$N_mz) = split(/\\:/,\$N_num_mz);	

###### correct mass before the search 		
				if(\$params->{mass_correction})
				{
					\$C_mz = \$C_mz - \$C_mz * \$mass_shift->{\$scan} / 1000000;
					\$N_mz = \$N_mz - \$N_mz * \$mass_shift->{\$scan} / 1000000;					
				}
				
					
				foreach my \$theoretical_formula (\@\$return) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
						if(\$N_num eq "N0")
						{
							if(\$theoretical_formula!~/N/)
							{
								\$nitrogen_match=1;
							}						
						}
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = 1;
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\\n";
						exit(1);
					}					
										
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{
						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
						\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
						\$formula_db=~s/(\\D+)1\$/\$1/g;						
						
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t", "pass", "\\t", \$formula_db, "\\t \\n";
						\$detected_formula++;						
					}
								
				}
				if(\$detected_formula==0 and \$params->{'adduct'} == 1)
				{
					my \$adduct_hash = \$p->get_adducts(\$params);
					foreach \$adduct_name (keys \%\$adduct_hash)
					{
						my \$mono_mass_no_adduct = \$mono_mass - \$adduct_hash->{\$adduct_name};
						my \$N_mz_no_adduct = \$C_mz_no_adduct = 0;
						if(\$N_mz>0)
						{
							\$N_mz_no_adduct = \$N_mz - \$adduct_hash->{\$adduct_name};
						}
						if(\$C_mz>0)
						{
							\$C_mz_no_adduct = \$C_mz - \$adduct_hash->{\$adduct_name};
						}					
						if(\$params->{'mode'} == -1)
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct + 1.007277;
						}
						else
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct - 1.007277;
						}					
						(\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass_no_adduct, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);	
						if(\$params->{'mode'} == -1)
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct - 1.007277;
						}
						else
						{
							\$mono_mass_no_adduct = \$mono_mass_no_adduct + 1.007277;
						}						
						
						foreach my \$theoretical_formula (\@\$return) 
						{
							chomp \$theoretical_formula;
							my \$carbon_match = 0;
							my \$nitrogen_match = 0;
							
							if(\$params->{'labeled_ID_method'} == 2)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
								if(\$N_num eq "N0")
								{
									if(\$theoretical_formula!~/N/)
									{
										\$nitrogen_match=1;
									}						
								}
							}
							elsif(\$params->{'labeled_ID_method'} == 3)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = 1;
							}
							elsif(\$params->{'labeled_ID_method'} == 4)
							{
								\$carbon_match = 1;
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\\\n";
								exit(1);
							}					
												
							if(\$carbon_match == 1 and \$nitrogen_match ==1)
							{
								
								my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
								my \$adduct_mass = \$adduct_hash->{\$adduct_name};
								print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t", "pass", "\\t", \$formula_db, "\\t",  \$adduct_name,"\\n";
								\$adduct_formula++;						
							}
										
						}
						
						foreach my \$theoretical_formula (\@\$return_norm) 
						{
						
						
							chomp \$theoretical_formula;				
							my \$carbon_match = 0;
							my \$nitrogen_match = 0;
							
							if(\$params->{'labeled_ID_method'} == 2)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
								if(\$N_num eq "N0")
								{
									if(\$theoretical_formula!~/N/)
									{
										\$nitrogen_match=1;
									}						
								}
							}
							elsif(\$params->{'labeled_ID_method'} == 3)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = 1;
							}
							elsif(\$params->{'labeled_ID_method'} == 4)
							{
								\$carbon_match = 1;
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\\n";
								exit(1);
							}					
								
							if(\$carbon_match == 1 and \$nitrogen_match ==1)
							{
							
								my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

								my \$query = new Spiders::DatabaseQuery();
								\$query->setBin("$lib");										
								\$query->setDatabase(\$params->{structure_database});
								\$query->setFormula(\$formula_db);
								\$query->setMass(\$mass_db);
								\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
								#\$query->setDatatype("SMILES");
								\$query->setDBname(\$params->{database});						
								\$str_db_return = \$query->QueryStructureDatabase();				
								my \$adduct = (\$adduct_formula > 0) ? \$adduct_name : "";
								my \$adduct_mass = (\$adduct_formula > 0) ? \$adduct_hash->{\$adduct_name} : 0;
								\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
								\$formula_db=~s/(\\D+)1\$/\$1/g;									
								print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Target","\\t",\$formula_db, "\\t",  \$adduct,"\\n";					
								
								if(scalar(\@\$str_db_return)>0)
								{
									foreach (\@\$str_db_return)
									{
										chomp \$_;							
										my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);

										print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  \$adduct,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
										
									}
								}
								
								\$index++;		
							}	
						}
			################################################

						
						
			###### generate decoys 
						foreach my \$theoretical_formula (\@\$return_decoy) 
						{
							chomp \$theoretical_formula;
							my \$carbon_match = 0;
							my \$nitrogen_match = 0;
							
							if(\$params->{'labeled_ID_method'} == 2)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
								if(\$N_num eq "N0")
								{
									if(\$theoretical_formula!~/N/)
									{
										\$nitrogen_match=1;
									}						
								}
							}
							elsif(\$params->{'labeled_ID_method'} == 3)
							{
								\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
								\$nitrogen_match = 1;
							}
							elsif(\$params->{'labeled_ID_method'} == 4)
							{
								\$carbon_match = 1;
								\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							}
							else
							{
								print "please set the labeled_ID_method parameter [1-4]\\n";
								exit(1);
							}					
								
							if(\$carbon_match == 1 and \$nitrogen_match ==1)
							{				
								
								my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

								my \$updated_formula_db = \$formula_db;
								if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
								{
									my \$Hminus = \$2 - 1  * \$decoy_H_num;
									\$updated_formula_db = \$1 . \$Hminus . \$3
								}

											
								my \$updated_mass_db = \$mass_db - 1.00782503207  * \$decoy_H_num;
									
								my \$query = new Spiders::DatabaseQuery();
								\$query->setBin("$lib");
								\$query->setDatabase(\$params->{structure_database});
								\$query->setFormula(\$updated_formula_db);
								\$query->setMass(\$updated_mass_db);
								\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
								#\$query->setDatatype("SMILES");						
								\$query->setDBname(\$params->{database});
								
								my \$str_db_return = \$query->QueryStructureDatabase();	
								my \$adduct = (\$adduct_formula > 0) ? \$adduct_name : "";
								\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
								\$formula_db=~s/(\\D+)1\$/\$1/g;									
								my \$adduct_mass = (\$adduct_formula > 0) ? \$adduct_hash->{\$adduct_name} : 0;                    						
								print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Decoy","\\t",\$formula_db, "\\t",  \$adduct,"\\n";
								
								if(scalar(\@\$str_db_return)>0)
								{
									foreach (\@\$str_db_return)
									{
										chomp \$_;
										my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);

										print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  \$adduct,"\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tDecoy","\\n";

										my \$smiles_orig = \$smiles;
										\$index++;

									}
								}						
								
								\$index++;		
							}	
						}
					}
	
				}
				else
				{

####### generating structure and scoring #######################
			
					foreach my \$theoretical_formula (\@\$return_norm) 
					{
					
					
						chomp \$theoretical_formula;				
						my \$carbon_match = 0;
						my \$nitrogen_match = 0;
						
						if(\$params->{'labeled_ID_method'} == 2)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							if(\$N_num eq "N0")
							{
								if(\$theoretical_formula!~/N/)
								{
									\$nitrogen_match=1;
								}						
							}
						}
						elsif(\$params->{'labeled_ID_method'} == 3)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = 1;
						}
						elsif(\$params->{'labeled_ID_method'} == 4)
						{
							\$carbon_match = 1;
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
						}
						else
						{
							print "please set the labeled_ID_method parameter [1-4]\\n";
							exit(1);
						}					
							
						if(\$carbon_match == 1 and \$nitrogen_match ==1)
						{
						
							my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

							my \$query = new Spiders::DatabaseQuery();
							\$query->setBin("$lib");										
							\$query->setDatabase(\$params->{structure_database});
							\$query->setFormula(\$formula_db);
							\$query->setMass(\$mass_db);
							\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
							#\$query->setDatatype("SMILES");
							\$query->setDBname(\$params->{database});						
							\$str_db_return = \$query->QueryStructureDatabase();				

							\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
							\$formula_db=~s/(\\D+)1\$/\$1/g;								
							print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Target","\\t",\$formula_db, "\\t",  \$adduct,"\\n";					
							
							if(scalar(\@\$str_db_return)>0)
							{
								foreach (\@\$str_db_return)
								{
									chomp \$_;							
									my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);

									print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  "","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tTarget","\\n";								
									
								}
							}
							
							\$index++;		
						}	
					}
		################################################

					
					
		###### generate decoys 
					foreach my \$theoretical_formula (\@\$return_decoy) 
					{
						chomp \$theoretical_formula;
						my \$carbon_match = 0;
						my \$nitrogen_match = 0;
						
						if(\$params->{'labeled_ID_method'} == 2)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
							if(\$N_num eq "N0")
							{
								if(\$theoretical_formula!~/N/)
								{
									\$nitrogen_match=1;
								}						
							}
						}
						elsif(\$params->{'labeled_ID_method'} == 3)
						{
							\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
							\$nitrogen_match = 1;
						}
						elsif(\$params->{'labeled_ID_method'} == 4)
						{
							\$carbon_match = 1;
							\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
						}
						else
						{
							print "please set the labeled_ID_method parameter [1-4]\\n";
							exit(1);
						}					
							
						if(\$carbon_match == 1 and \$nitrogen_match ==1)
						{				
							
							my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

							my \$updated_formula_db = \$formula_db;
							if(\$formula_db=~/(\.\*H)(\\d\+)(\\D\+\.\*)/)
							{
								my \$Hminus = \$2 - 1  * \$decoy_H_num;
								\$updated_formula_db = \$1 . \$Hminus . \$3
							}

										
							my \$updated_mass_db = \$mass_db - 1.00782503207  * \$decoy_H_num;

							my \$query = new Spiders::DatabaseQuery();
							\$query->setBin("$lib");
							\$query->setDatabase(\$params->{structure_database});
							\$query->setFormula(\$updated_formula_db);
							\$query->setMass(\$updated_mass_db);
							\$query->setDatatype("INCHIKEY,SMILES,GENERALNAME");
							#\$query->setDatatype("SMILES");						
							\$query->setDBname(\$params->{database});
							
							my \$str_db_return = \$query->QueryStructureDatabase();	
							\$formula_db=~s/(\\D+)1(\\D+)/\$1\$2/g;
							\$formula_db=~s/(\\D+)1\$/\$1/g;	                  						
							print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t","Decoy","\\t",\$formula_db, "\\t",  \$adduct,"\\n";
							
							if(scalar(\@\$str_db_return)>0)
							{
								foreach (\@\$str_db_return)
								{
									chomp \$_;
									my (\$formula_db,\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);

									print RESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_comb->{\$N_mz}->{'intensity'},"\\t",\$formula_comb->{\$C_mz}->{'intensity'},"\\t",\$pair_score,"\\t",\$formula_db,"\\t",  "","\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\tDecoy","\\n";

									my \$smiles_orig = \$smiles;
									\$index++;

								}
							}						
							
							\$index++;		
						}	
					}
				}
			}			
		}	
    }
	
}
		
EOF
	close(RUNSHELL);	
}




sub Check_Job_stat
{
	my ($self,$jobs_prefix,$job_num,$dta_path) = @_;
	my $params = $self->get_parameter();	
	my $LOG = $self->get_log_file();		
	my $job_info=1;
    my ($username) = getpwuid($<);
	my $command_line="";
	my $dot = ".";
	while($job_info)
	{
		if($params->{'job_management_system'} eq 'LSF')
		{
			$command_line =  "bjobs -u $username";
		}
		elsif($params->{'job_management_system'} eq 'SGE')
		{
			$command_line =  "qstat -u $username";
		}


		#Consider only the one that we submitted
		if($params->{'job_management_system'} eq 'LSF')
		{
			$command_line =  "bjobs -u $username";	

			my $job_status=qx[$command_line];
			my @job_status_array=split(/\n/,$job_status);
			my $job_number = $job_num - scalar (@job_status_array) + 1;
			if(scalar (@job_status_array) == 0)
			{
				print "\r  $job_num jobs finished          ";
			}
			else
			{
				print "\r  $job_number jobs finished          ";
				sleep(5);
			}
			if(scalar(@job_status_array)>0)
			{
				$job_info=1;				
			}
			else
			{
				$job_info=0;		
			}			
		}
		elsif($params->{'job_management_system'} eq 'SGE')
		{
			my $job_status=qx{$command_line 2>&1};
			
			my @job_status_array=split(/\n/,$job_status);		
			@job_status_array = grep(/$jobs_prefix/,@job_status_array);
			if($jobs_prefix =~ /pair_/)
			{			
				if($job_status=~/No unfinished job found/)
				{
					$job_info=0;
					print "  \n";
				}
				elsif((scalar(@job_status_array))==0)
				{
					$job_info=0;
				}
				else
				{
					sleep(30);
				}
			}	
		
			elsif($jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/)
			{
				my $check_command = "ls -f $dta_path\/\*.smout \| wc -l";
				my @outfile = glob("$dta_path\/\*.smout");

				my $outfile_num=scalar @outfile;

				print "\r  $outfile_num files have done         ";
				sleep(30);
				if((scalar(@job_status_array))==0)
				{
					$job_info=0;
				}				
	
			}
			elsif($jobs_prefix =~ /frag_/)
			{
				if($params->{'job_management_system'} eq 'LSF')
				{	
					$command_line =  "bjobs -u $username";	
				}
				elsif($params->{'job_management_system'} eq 'SGE')
				{
					$command_line =  "qstat -u $username";			
				}			
				my $job_status=qx[$command_line];
				my @job_status_array=split(/\n/,$job_status);
				my $job_number = $job_num - scalar (@job_status_array) + 2;
				if(scalar (@job_status_array) == 0)
				{
					print "\r  $job_num jobs finished          ";					
				}
				else
				{
					print "\r  $job_number jobs finished          ";
					sleep(30);
				}
				if(scalar(@job_status_array)>0)
				{
					$job_info=1;				
				}
				else
				{
					$job_info=0;		
				}
			}
			
		}

	}
}


sub LuchParallelJob{
	
	my($self,$FileName,$cmd,$GridType,$outputName,$dta_path)= @_;
	
	open(JOB,">$FileName") || die "can not open $FileName\n";
	if($GridType eq 'LSF')
	{	
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q standard\n";
			print JOB "#BSUB -M 20000\n";
			print JOB "#BSUB -R \"rusage[mem=10000]\"\n";			
			
			print JOB "#BSUB -eo $dta_path/$outputName.e\n";
			print JOB "#BSUB -oo $dta_path/$outputName.o\n";
			print JOB $cmd;		
			close(JOB);
			system(qq(bsub <$FileName >/dev/null 2>&1));	
	}
	if($GridType eq 'SGE')
	{

		print JOB "#!/bin/bash\n";
		print JOB "#\$ \-S /bin/bash\n";  #In our cluster this line is esential for executing some bash commands such as for
		print JOB "#\$ \-N $outputName\n";
		print JOB "#\$ \-e $dta_path/$outputName.e\n";
		print JOB "#\$ \-o $dta_path/$outputName.o\n";
		print JOB $cmd;
		close(JOB);
		system(qq(qsub -cwd -cwd -pe mpi 4 -l mem_free=8G,h_vmem=6G $FileName >/dev/null 2>&1));
	}
	if($GridType eq 'PBS')
	{

		print JOB "#!/bin/bash\n";
		print JOB "#PBS -N $outputName\n";
		print JOB "#PBS -e $dta_path/$outputName.e\n"; 
		print JOB "#PBS -o $dta_path/$outputName.o\n"; 			
		print JOB $cmd;
		close(JOB);
		system(qq(qsub -cwd $FileName >/dev/null 2>&1));
	}
	close(JOB);
}

sub submit_jobs
{
	my ($self,$job_num,$job_name,$working_path) = @_;

	my $params = $self->get_parameter();
	my $LOG = $self->get_log_file();
	#my $MAX_PROCESSES = $params->{'processor_number'};	
	my $MAX_PROCESSES = 4;
	######### running jobs ######################### 
	my $job_list;
	if($params->{'fragment_method'}==2 and $job_name=~/frag/)
	{
		$params->{'cluster'} = 0;
	}
	if($params->{'cluster'} eq '1')
	{
		if($params->{'job_management_system'} eq 'LSF')
		{
			for(my $i=0;$i<$job_num;$i++)
			{
				my $command_line = qq(cd $working_path && bsub <${job_name}_${i}.sh);
				my $job=qx[$command_line];
				chomp $job;
				my $job_id=0;
				if($job=~/Job \<(\d*)\> is/)
				{
					$job_id=$1;
				}
				$job_list->{$job_id}=1;
			}
		}
		elsif($params->{'job_management_system'} eq 'SGE')
		{

			for(my $i=0;$i<$job_num;$i++)
			{
				my $job_name = "${job_name}_${i}.sh";
				my $command_line = qq(cd $working_path && qsub -cwd -pe mpi 2 $job_name);
				my $job=qx[$command_line];
				chomp $job;
				my $job_id=0;
				if($job=~/$job_name \<(\d*)\> is/)
				{
					$job_id=$1;
				}
				$job_list->{$job_id}=1;
				my $count = $i+1;
				print "\r  $count jobs were submitted";				
			}
			print  $LOG "  $job_num jobs were submitted\n";	
		}
		print "\n";
#		print "\n  You submitted $job_num jobs for database search\n";
#		print LOG "\n  You submitted $job_num jobs for database search\n";		
		$self->Check_Job_stat("${job_name}_",$job_num,$working_path);		
	}
	elsif($params->{'cluster'} eq '0')
	{
		print "  It may take several hours to run jobs using a single server. please be patient!\n";
		print $LOG "  It may take several hours to run jobs using a single server. please be patient!\n";		
        my $pm = new Parallel::ForkManager($MAX_PROCESSES);
		my $job_finish=0;
#		$pm->run_on_wait(sub{$job_finish++;print "\r  $job_finish jobs finished";},0.5);
        for my $i ( 0 .. $job_num )
        {
            $pm->start and next;
			my $job_name = "${job_name}_${i}.sh";			
            system("cd $working_path && sh $job_name >/dev/null 2>&1");
            $pm->finish; # Terminates the child process
#			print "\r  $i jobs finished";				
        }
		#Check_Job_stat("${job_name}_",$job_num,$dta_path);
	#	print  $LOG "  $MAX_PROCESSES jobs were submitted\n";					
        $pm->wait_all_children;		
	}

}


sub runjobs
{
	my ($self,$file_array,$dta_path,$job_name,$ratio,$defect,$parameter) = @_;
	my $params = $self->get_parameter();	
	my $curr_dir = getcwd();
    	#my $MAX_PROCESSES = $params->{'processor_number'};
    	my $MAX_PROCESSES = 4;	
	my $job_num = 400;
	my $shell = "";
	if($job_name =~ /pair/)
	{
		$shell = "pair_shell.pl"
	}
	elsif($job_name =~ /sch/)
	{
		$shell = "runsearch_shell.pl"
	}

	
	my $dta_num_per_file = 10;
	if($params->{'cluster'} eq '0')
	{
		$job_num = $MAX_PROCESSES;
		$dta_num_per_file = int($#$file_array / $job_num) + 1;	
	}
	else
	{
		$job_num = int($#$file_array / $dta_num_per_file) + 1;
	}
	## Set the maximum number of jobs to 400
	if ($job_num > 400) {
		$job_num = 400;
		$dta_num_per_file = int($#$file_array / $job_num) + 1;
	}	
	

 
	for(my $i = 0; $i < $job_num; $i++)
	{	
		if (($i * $dta_num_per_file) > $#$file_array) 
		{
			$job_num = $i;
			last;
		}	
		
		open(JOB,">$dta_path/${job_name}_${i}.sh") || die "can not open the job files\n";

		my @dta_file_arrays=();
		my $multiple_jobs_num  = 0;
		for(my $j=0;$j<$dta_num_per_file;$j++)
		{
			if(($i*$dta_num_per_file+$j)<=$#$file_array)
			{			
				push (@dta_file_arrays,$$file_array[$i*$dta_num_per_file+$j]);
			}
		}
				
		if($params->{'job_management_system'} eq 'LSF')
		{
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q standard\n";
			print JOB "#BSUB -R \"rusage[mem=10000]\"\n";
			print JOB "#BSUB -eo $dta_path/${job_name}_${i}.e\n";
			print JOB "#BSUB -oo $dta_path/${job_name}_${i}.o\n";
	
		}
		elsif($params->{'job_management_system'} eq 'SGE')
		{
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -N ${job_name}_${i}\n";
			print JOB "#\$ -e ${job_name}_${i}.e\n";
			print JOB "#\$ -o ${job_name}_${i}.o\n";			
		}
		if($job_name =~ /pair/)
		{
			foreach (@dta_file_arrays)
			{
				next if(!defined($_));
				print JOB "perl $dta_path/$shell -param $parameter $_ \n";	
			}			
		}
		elsif($job_name =~ /frag/)
		{
			foreach (@dta_file_arrays)
			{
				next if(!defined($_));
				print JOB "$_ \n";	
			}			
		}			
		elsif($job_name =~ /sch/)
		{
		
			if($params->{'labeled_data'} eq '0')
			{
				foreach (@dta_file_arrays)
				{
					next if(!defined($_));			
					print JOB "perl $dta_path/$shell -param $parameter $_ \n";
				}					
			}
			else
			{
				foreach (@dta_file_arrays)
				{
					if($defect eq "Undef")
					{
						next if(!defined($_));
						print JOB "perl $dta_path/$shell -param $parameter $_ \n";						
					}
					elsif(defined $defect)
					{
						next if(!defined($_));			
						print JOB "perl $dta_path/$shell -param $parameter $_ -NC $ratio->{'NC'} -CC $ratio->{'CC'} -NC_std $ratio->{'NC_std'} -CC_std $ratio->{'CC_std'} -NC_defect_loc $defect->{'NC'}->[0] -CC_defect_loc $defect->{'CC'}->[0] -NC_defect_scale $defect->{'NC'}->[1] -CC_defect_scale $defect->{'CC'}->[1]\n";	
					}
					else
					{
						next if(!defined($_));
						print JOB "perl $dta_path/$shell -param $parameter $_ \n";						
					}
				}
			
			}			
			
		}		
		
		close(JOB);
	}
	my $total_file_num = $#$file_array + 1;
	print "  Generating $job_num jobs containing $total_file_num files\n";
	$self->submit_jobs($job_num,$job_name,$dta_path);
}

sub luanch_ms2_jobs
{
	my ($self,$msms_hash_ref,$ms_hash_mol,$dta_path,$parameter) = @_;

	my $count = 0;
	my $MS1_MS2_matched;
	my @frag_job_list=();
	my %msms_hash = %$msms_hash_ref;
	my $params = $self->get_parameter();
	my $lib = $self->get_library_path();
	my %ms2_scan_num_hash;
	foreach my $scan (sort {$a<=>$b} keys %msms_hash)
	{

		print "\r  processing MS2 scan: $scan";
		my $ms_prec_mh = $msms_hash{$scan}{'prec_MH'};

		next if(!defined($ms_prec_mh));		
		foreach my $scan_missile (keys %$ms_hash_mol)
		{

			next if($scan < ($scan_missile-$params->{'matched_scan_dist'} ) or $scan > ($scan_missile+$params->{'matched_scan_dist'}));
			next if(!defined($ms_hash_mol->{$scan_missile}));
			#my $i=0;
			foreach my $mz_missile (keys %{$ms_hash_mol->{$scan_missile}})
			{

				my $tolerance_Da = $params->{'isolation_window'};				
				if(abs($mz_missile-$ms_prec_mh)<$tolerance_Da)
				{
					$ms2_scan_num_hash{$scan} = 1;
					foreach my $smile (keys %{$ms_hash_mol->{$scan_missile}->{$mz_missile}})
					{
						$MS1_MS2_matched->{$scan_missile}->{$mz_missile}->{$scan} = 1; 
						my ($prec_type) = keys %{$ms_hash_mol->{$scan_missile}->{$mz_missile}->{$smile}};
						my $lines = $ms_hash_mol->{$scan_missile}->{$mz_missile}->{$smile}->{$prec_type};
						chomp $lines;
						my @data=split(/\t/,$lines);
						
						my $target_decoy = $data[14];
						$smile =~s/\(/\\\(/g;
						$smile =~s/\)/\\\)/g;
						$smile ="\"" . $smile . "\"";

						#print "\n",$scan,"\t",$scan_missile,"\t",$smile,"\t",$mz_missile,"\t",$ms_prec_mh,"\t",$scan_missile,"\n";
						#push (@frag_job_list,"perl $lib/frag_shell.pl -dtafile $dta_path/${scan}.MS2  -smile $smile -mass 50 -depth 2 -param  $parameter -ptype $prec_type -prec_mh $ms_prec_mh -MS1_scan $scan_missile -ms1_mass $mz_missile");
						push (@frag_job_list,"perl $lib/frag_shell.pl -dtafile $dta_path/${scan}.MS2  -smile $smile -mass 50 -depth 2 -param  $parameter -ptype $prec_type -target_decoy $target_decoy");
		#				$i++;
																								
					}
				}
				#last if($i>1);							
			}
			#last if($i>1);
		}		
	}
	my $ms2_scan_num = scalar keys (%ms2_scan_num_hash);
	return (\@frag_job_list,$MS1_MS2_matched,$ms2_scan_num);
}



1;
