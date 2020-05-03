--------------------- 
Contents of this file 
--------------------- 
 * Introduction 
 * Software Requirements 
 * Hardware Requirements 
 * Installation 
 * Database Setup 
 * Command Line Arguments 
 * Testing JUMPm Using Example Files 
 * Output 
 * Contact 
 
------------ 
Introduction 
------------ 
JUMPm is a program for untargeted metabolite identification for liquid chromatography and tandem mass spectrometry metabolomics. This program determines chemical formulas from either unlabeled or stable-isotope labeled metabolomics data, and derives possible structures by predictive fragmentation during database search. JUMPm uses a target-decoy strategy based on the octet rule to estimate the rate of false discovery (FDR).  The user specifies a target FDR and JUMPm will filter the data to reach the target.  FDR is a critical measure of confidence which researchers can use in the analysis of their data. The program is written in perl and is designed for high performance parallel computing systems.
  
---------------------
Software Requirements
---------------------
The program is written in a combination of Perl, Java and R. It should run on every system with a Perl5, Java 1.7 (or higher) interpreter and R 3.1.0 (or higher).
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


R package required (for unlabeled data):
        mzR (biocoductor)

To install mzR package, please visit http://biodonductor.org
 
--------------------- 
Hardware Requirements 
---------------------  
The program can be run on either high performance computing systme or a single server.  
  
To run on a cluster: 
  Batch-queuing system: SGE is recommended (version 6.1u5, qsub and qstat) 
  32 GB memory on each node 
 
To run on a single server 
  32 GB memory 
  2 GHz CPU processors with a minimum of 4 cores 
   
------------
Installation
------------
The source code, database, sample data can be download from https://www.stjuderesearch.org/site/lab/peng/jumpm and
https://github.com/JUMPSuite/JUMPm (source code only). Database includes both mass formula and structure databases.
After downloading the source code, you can put it in any working directory (e.g. /home/xxxx/JUMPm).
*If you install the program on the cluster, the folder containing all source codes is required to be accessible from each node.
For example, you can put it on your local home directory.

--------------
Database Setup
--------------
The database can be downloaded from http://www.stjuderesearch.org/site/lab/peng, including mass formula database and structure database.
After downloading databases, please uncompress the two databases on your local directory and set the two paths in the parameter file.

MASS_FORMULA_DB: http://ftp.stjude.org/pub/software/JUMPm/MASS_FORMULA_DB.tar.gz
STRUCTURE_DB: http://ftp.stjude.org/pub/software/JUMPm/STRUCTURE_DB.tar.gz

----------------------
Command Line Arguments
----------------------
Example:  perl <JUMPm Installation Path>/jumpm -p <JUMPm parameter file> <mzXML file(s)>
1. Please use the absolute path for jumpm program. For example: perl /home/user/JUMPm/jumpm
2. -p <file> specifies a JUMP parameter file; please keep the parameter in the same folder of the mzXML file(s) and use relative path instead of absolute path
3. <mzXML file(s)> specifies one or multiple mzXML file(s); please use relative path instead of absolute path

An example parameter file is provided in the package of source code.
JUMPm assumes mzXML file(s) is/are converted from .RAW file(s) using ReAdW (without compression, 32-bit precision encoding/decoding).
 
---------------------------------
Testing JUMPm using example files
---------------------------------
1. Two LC-MS/MS files in .mzXML format are provided in JUMPm package with folder name of Example
        Example
        |-----labeled_data
                |------- yeast_4plex_HILIC_NEG_ph3_HM.mzXML
                |------- jumpm_labeled_test.params
        |-----unlabeled_data
                |------ 12C_HILIC_Neg_1.mzXML
                |------ jumpm_unlabeled_test.params

To test JUMPm program, we only limit the scan range 1000 to 1050, which is defined in the parameter files.
If you want to search entire run, please change the following two parameters to cover all scans.
        first_scan_extraction = 0
        last_scan_extraction = 10000000

2. If you want to set up JUMPm program on a single server Linux system, please download the JUMPm package by the following command:
        2.1 move to the directory where you want to install JUMPm (e.g. cd /home/usr/jumpm)
        2.2 wget http://ftp.stjude.org/pub/software/JUMPm/JUMPm.tar.gz
        2.3 uncompress the source code: tar -vxzf JUMPm.tar.gz
        2.4 if all modules/packages associated with the program are installed appropriately, please test the program
 
------
Output
------
For labeled dataset(s),
1. metabolite feature
This table lists all features detected in each LC-MS run, which only contains information about features before performing database search. All features listed in this table are de-isotoped features. In other words, isotopic features are excluded in this table.
The format of metabolite feature table is shown below:
index   m/z     z       MS1 scan#       Intensity       S/N

2. formula table
This table lists all formulas detected for each each feature. The table lists the pairing information for a labeled dataset, whereas it use "N/A" for three columns, including N15 mass, C13 mass and Pscore. The target/decoy column shows whether the formula is a target or decoy hits. Two levels of target/decoy are used: the first type is only applied to the hits from user-specified database; the second is applied to all hits from theorectical formula database based on the octet rule.
The format of metabolite feature table is shown below:
index   MS1     C12 mass        N15 mass        C13 mass        C12 intensity   N15 intensity   C13 intensity   S/N     Pscore  Formula Target/Decoy

3. Spectrum-match table
This table lists all structure candidates for each formula. Three structure related features are added to this table: smiles, chemical name and Mscore.
The format of metabolite feature table is shown below:
index   MS1     C12 mass        N15 mass        C13 mass        C12 intensity   N15 intensity   C13 intensity   S/N     Pscore  Formula Target/Decoy    structure (SMILES)      name    Mscore

4. Structure table
This table lists a unique structure candidate with the best Mscore for each formula.
The format of metabolite feature table is shown below:
index   MS1     C12 mass        N15 mass        C13 mass        C12 intensity   N15 intensity   C13 intensity   S/N     Pscore  Formula Target/Decoy    structure (SMILES)      name    Mscore

5. Intermediate files
All intermediate files can be found in the intermediate folder, called RAW_FILE_NAME.1.intermediate and RAW_FILE_NAME.1.

For unlabeled dataset(s),
1. result file
This table lists all compounds identified by the features from run(s). Not to lose any information, all possible mappings between candidate compounds and features are listed with Mscore. The format of the table is as follows:
(nominal)FeatureNo      Ion     m/z     RT      Formula Target/Decoy    Name    SMILES  InChiKey        Mscore  intensity in each run

-------
Contact
-------
* To submit bug reports and feature suggestions, or to track changes, please contact:
Xusheng Wang (xusheng.wang@und.edu), Ji-Hoon Cho (ji-hoon.cho@stjude.org) or Junmin Peng (junmin.peng@stjude.org)
