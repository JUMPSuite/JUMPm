----------------------
Contents of this file
----------------------

 * Introduction
 * Software Requirements
 * Hardware Requirements
 * Installation
 * Database Setup
 * Command Line Arguments
 * Testing JUMPm Using Example Files
 * Output
 * Maintainers

----------------------
Introduction
----------------------

JUMPm is a program for untargeted metabolite identification for liquid chromatography and tandem mass spectrometry metabolomics. The computer algorithm determines chemical formulas from either unlabeled or stable-isotope labeled metabolomics data, and derives possible structures by predictive fragmentation during database search. JUMPm uses a target-decoy strategy based on the octet rule to estimate the rate of false discovery (FDR).  The user specifies a target FDR and JUMPm will filter the data to reach the target.  FDR is a critical measure of confidence which researchers can use in the analysis of their data. The program is written in perl and is designed for high performance parallel computing systems.
 
----------------------
Software Requirements
---------------------- 

The program is written in a combination of Perl, Java and R. It should run on every system with a Perl5, Java 1.7 interpreter and R 3.1.0. 
The minimum required Perl version should be Perl 5.8 or better.

Below are all perl modules needed. Please make sure you have install all of these modules.

    Parallel::ForkManager
	Class::Std
	Statistics::R
	Statistics::Basic
	Statistics::Descriptive
	Set::Partition
	Regexp::Common
	Number::Format

To install perl modules, you can use cpan. please refer to:
http://www.cpan.org/modules/INSTALL.html
 

Java module: 
CDK

1. Make sure that you have installed JAVA with the latest version
2. CDK module is not required to be installed. 
	

The R source code (R-3.1.0.tar.gz) is provided in the JUMPm folder. You need to rebuild it as follows:

1. tar -vxf R-3.1.0.tar.gz
2. ./configure --with-readline=no
3. make

Note: You do not need to install the R because JUMPm called the R from JUMPm folder.
if you use default R environment, please change the $R_library variable in the main JUMPm script. 
 
----------------------
Hardware Requirements
---------------------- 

The program can be run on either high performance computing systme or a single server. 
 
To run on a cluster:
 Batch-queuing system: SGE, version 6.1u5, qsub and qstat
 32 GB memory on each node

To run on a single server
  32 GB memory
  2 GHz CPU processors with a minimum of 4 cores
  
----------------------
Installation
---------------------- 

The source code, database, sample data can be download from http://www.stjuderesearch.org/site/lab/peng/jumpm

Database includes MASS_FORMULA_DB and STRUCTURE_DB.

After downloading the source code, you can put it in any working directory (e.g. /home/xxxx/JUMPm). 
**If you install the program on the cluster, the folder containing all source code is required to be accessible to each node. For example, you can put it on your local home directory. 

----------------------
Database Setup
---------------------- 

The database can be downloaded from http://www.stjuderesearch.org/site/lab/peng, including mass formula database and structure database.
After downloading databases, you can uncompress the two databases on your local directory and set the two paths in the parameter file.  

MASS_FORMULA_DB: http://ftp.stjude.org/pub/software/JUMPm/MASS_FORMULA_DB.tar.gz
STRUCTURE_DB: http://ftp.stjude.org/pub/software/JUMPm/STRUCTURE_DB.tar.gz

----------------------
Command Line Arguments
----------------------

Example:  perl <JUMPm Installation Path>/jumpm -p <JUMPm parameter file> <MS/MS data file(s)>
Please use the absolute path for jumpm program. For example: perl /home/user/JUMPm/jumpm 
-p <file> specifies a JUMP parameter file; please keep the parameter in the same folder of the raw file and use relative path instead of absolute path
<MS/MS data file(s)> specifies one or multiple MS/MS data file(s); please use relative path instead of absolute path

An example parameter file is provided in the package of source code. 


MS/MS data file(s) can be either .RAW or mzXML file. 
The mzXML file is converted from .RAW file by ReAdW. A latest version (4.3.1) is included in the JUMPm package. 

To run the ReAdW on windows system, please follow the steps below: 

1. open a Command Prompt window using Start menu->run, enter cmd

2. change directory to ReAdW folder

3. enter the following command:
	regsvr32 XRawfile2.dll
	
4. enter the following command to convert RAW into mzXML

	ReAdW --mzXML C:\test\input.raw c:\test\output.mzXML

5. upload output.mzXML to a server with Linux system

----------------------
Testing JUMPm using example files
----------------------

1. Two LC-MS/MS files in .mzXML format are provided in JUMPm package with folder name of Example
	Example
	|-----labeled_data
	       |------- yeast_4plex_HILIC_NEG_ph3_HM.mzXML
		   |------- jumpm_negative_single_server.params
		   |------- jumpm_negative_cluster.params
	|-----unlabeled_data
			|------ 12C_HILIC_Neg_1.mzXML
			|------ jumpm_negative_single_server.params
			|------ jumpm_negative_cluster.params
			
To test JUMPm program, we only limit the scan range 1000 to 1050, which is defined in the parameter files.
If you want to search entire run, please change the following two parameters to cover all scans.
first_scan_extraction = 0 
last_scan_extraction = 10000000

For example, to test the labeled data and if you JUMPm program is installed in a folder of "/home/userA/JUMPm_v1.8.0/", please to the folder "/home/userA/JUMPm_v1.8.0/Example/labeled_data/". If you use a single server, 

perl /home/userA/JUMPm_v1.8.0/jumpm -p jumpm_negative_single_server.params yeast_4plex_HILIC_NEG_ph3_HM.mzXML

it may takes 1 hour to finish the run depending on your hardware and the range of scans. To boost the processing time, please reduce the scan range, for example
first_scan_extraction = 500 
last_scan_extraction = 800


2. If you want to set up JUMPm program on a single server Linux system, please download the JUMPm package by the following command:
	2.1 wget http://ftp.stjude.org/pub/software/JUMPm/JUMPm_v1.8.0.tar.gz
	2.2 uncompress the source code: tar -vxzf JUMPm_v1.0.2.tar.gz
	2.3 single server: if all modules associated with program are installed appropriately, you can go the the example folder (e.g. /home/user/jumpm_foler/labeled_data) and execute the program
		perl /home/user/jumpm_foler/jumpm -p jumpm_negative_single_server.params Yeast_4Plex_HILIC_Neg_1.mzXML
	2.4 cluster: please run the program as follows:
		perl /home/user/jumpm_foler/jumpm -p jumpm_negative_cluster.params Yeast_4Plex_HILIC_Neg_1.mzXML		


----------------------
Output
----------------------
1. metabolite feature
This table lists all features detected in each LC-MS run, which only contains information about features before performing database search. All features listed in this table are de-isotoped features. In other words, isotopic features are excluded in this table.
The format of metabolite feature table is shown below:
index	m/z	z	MS1 scan#	Intensity	S/N


2. formula table
This table lists all formulas detected for each each feature. The table lists the pairing information for a labeled dataset, whereas it use "N/A" for three columns, including N15 mass, C13 mass and Pscore. The target/decoy column shows whether the formula is a target or decoy hits. Two levels of target/decoy are used: the first type is only applied to the hits from user-specified database; the second is applied to all hits from theorectical formula database based on the octet rule.   
The format of metabolite feature table is shown below:
index	MS1	C12 mass	N15 mass	C13 mass	C12 intensity	N15 intensity	C13 intensity	S/N	Pscore	Formula	Target/Decoy

3. Spectrum-match table
This table lists all structure candidates for each formula. Three structure related features are added to this table: smiles, chemical name and Mscore.
The format of metabolite feature table is shown below:
index	MS1	C12 mass	N15 mass	C13 mass	C12 intensity	N15 intensity	C13 intensity	S/N	Pscore	Formula	Target/Decoy	structure (SMILES)	name	Mscore

4. Structure table
This table lists a unique structure candidate with the best Mscore for each formula. 
The format of metabolite feature table is shown below:
index	MS1	C12 mass	N15 mass	C13 mass	C12 intensity	N15 intensity	C13 intensity	S/N	Pscore	Formula	Target/Decoy	structure (SMILES)	name	Mscore

5. Intermediate files
All intermediate files can be found in the intermediate folder, called RAW_FILE_NAME.1.intermediate and RAW_FILE_NAME.1.


----------------------
Maintainers
----------------------

* To submit bug reports and feature suggestions, or to track changes, please contact:

 Xusheng Wang (xusheng.wang@stjude.org) and Junmin Peng (junmin.peng@stjude.org)
 
 