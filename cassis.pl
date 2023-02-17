#!/usr/bin/perl -w

###############################################################################
# Copyright (C) 2010 Université Claude Bernard Lyon 1                         #
#                                                                             #
# Contributors : Christian Baudet, Claire Lemaitre                            #
#                                                                             #
# Contact : christian.baudet@gmail.com, claire.lemaitre@gmail.com             #
#                                                                             #
# This file is part of Cassis.                                                #
#                                                                             #
# Cassis is free software: you can redistribute it and/or modify              #
# it under the terms of the GNU General Public License as published by        #
#  the Free Software Foundation, either version 3 of the License, or          #
# (at your option) any later version.                                         #
#                                                                             #
# Cassis is distributed in the hope that it will be useful,                   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with Cassis.  If not, see <http://www.gnu.org/licenses/>.             #
###############################################################################

use strict;
use File::Spec;
use FindBin qw($Bin);
require "${Bin}/core/util.pm";

###############################################################################
## CONSTANTS ##################################################################
###############################################################################

# GENERAL PARAMETERS - DEFAULT VALUES
my $MAXIMUM_SEQUENCE_R_SIZE     = 1000000000;
my $MAXIMUM_SEQUENCE_AB_SIZE    = 3000000;
my $EXTEND_BEFORE_VERIFY_LENGTH = "T";

# GENES - DEFAULT VALUES
my $EXTEND_GENES               = "T";
my $EXTEND_AB_GENES            = "T";
my $MINIMUM_SEQUENCE_SIZE_GENE = 1;

# BLOCKS - DEFAULT VALUES
my $EXTEND_BLOCKS                = 50000;
my $EXTEND_AB_BLOCKS             = 50000;
my $MINIMUM_SEQUENCE_SIZE_BLOCKS = 50000;

# LASTZ - DEFAULT VALUES
my $LASTZ_LEVEL             = 2;
my $LASTZ_PATH              = "${Bin}/lastz/lastz";
my $LASTZ_LEVEL1_PARAMETERS = "K=3000 E=30 O=400 L=2200 H=2000 B=2";
my $LASTZ_LEVEL2_PARAMETERS = "K=3000 E=30 O=400 L=3000 H=2000 B=2";
my $LASTZ_LEVEL3_PARAMETERS = "K=2200 L=6000 H=2000 B=2";
my $LASTZ_LEVEL1_MATRIX     = "${Bin}/lastz/matrix-level1.txt";
my $LASTZ_LEVEL2_MATRIX     = "default";
my $LASTZ_LEVEL3_MATRIX     = "${Bin}/lastz/matrix-level3.txt";

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################

# Verify the environment
my %ENVIRONMENT;
verifyEnvironment();

# Receive and process the parameters
my %PARAMETERS;
validateParameters();

# 1st step: create the directory structure
createWorkDirectories();

# 2nd step: validate the FASTA directories
validateFastaDirectories();

# 3rd step: validate table of genes / synteny blocks
validateInputTable();

# 4th step: identify the breakpoint regions
identifyBreakpoints();

# 5th step: create the FASTA files (sequences SR, SA and, SB)
createBreakpointFASTAs();

# 6th step: align the sequences
alignBreakpoints();

# 7th step: create dotplot
dotplotBreakpoints();

# 8th step: perform segmentation
refineBreakpoints();

exit(0);

###############################################################################
## AUXILIARY SUBROUTINES ######################################################
###############################################################################

###############################################################################
# This function prints the script usage
sub printUsage {

	my ($error) = @_;

	if ( defined $error ) {
		print STDERR "\nERROR: $error\n\n";
	}

	print << "End_print_usage";

Usage:
  perl cassis.pl [options] <table> <type> <dirGR> <dirGO> <outputdir>

-------------------------------------------------------------------------------
Mandatory Parameters:

  <table> Path for the file that contains the table of orthologous genes or
          synteny blocks.

  <type> Type of the table: G for orthologous genes and B for synteny blocks.

  <dirGR> Path for the directory where the script can locate the FASTA files 
          of the chromosomes of the reference genome (GR).

  <dirGO> Path for the directory where the script can locate the FASTA files 
          of the chromosomes of the genome (GO) which will be compared with  
          the reference genome.

  <outputdir> Name of the directory where the script will write the results.
              The directory must exist (the script will not try to create it).

-------------------------------------------------------------------------------
Optional parameter:

  --lastzlevel N

  LASTZ level: N = 1, 2 or 3
  + Level 1 : LASTZ parameters for closely related species.
  + Level 2 : Default level
  + Level 3 : LASTZ parameters for distantly related species.
  [Default value = 2]


  --max_length_sr N

  Maximum length for sequence SR: N = integer bigger than zero
  If SR length is bigger than N, the breakpoint receives the
  status -6 and it is not processed by the segmentation step.
  [Default value = 1000000000]


  --max_length_sab N

  Maximum length for sequences SA and SB: N = integer bigger than zero
  If SA (or SB) length is bigger than N, the length is corrected
  to fit this limit value.
  [Default value = 3000000]
 

  --minimum_length N

  Minimum sequence length: N = integer bigger than zero
  If one or more sequences have length smaller than N, the breakpoint receives
  one of the status -2, -3 , -4 or -5 and  it is not processed by the 
  segmentation step.
  [Default value = 1 (for genes) or 50000 (for synteny blocks)]


  --extend_before [T/F]

  Extend before verify length: True [T] or FALSE [F]
  This parameter determines if the method verifies the minimum sequence length
  before (F) or after (T) extending the sequence.
  [Default value = T]


  --extend_by_adding_gene [T/F]

  Extend sequences SR, SA and SB: True [T] or FALSE [F]
  Extend sequences SR, SA and SB by adding the orthologous genes which are in
  the boundaries of them. (Orthologous genes: pairs (Ar,Ao) and (Br,Bo))
  Warning: This parameter is available only for table of orthologous genes.
  [Default value = T]


  --extend_ab_by_adding_gene [T/F]

  Extend sequences SA and SB: True [T] or FALSE [F]
  Extend sequences SA and SB by adding the non orthologous genes which are in
  the boundaries of them. (Non orthologous genes: Co and Do)
  Warning: This parameter is available only for table of orthologous genes.
  [Default value = T]


  --extend_by_adding_fragment N

  Extend sequences SR, SA and SB: N = integer bigger than or equal to zero
  Extend sequences SR, SA and SB by adding to them a fragment of length N 
  from the sides that are supposed to be orthologous. It is equivalent to add
  the orthologous genes (Ar,Ao) and (Br,Bo) to the sequences.
  Warning: This parameter is available only for table of synteny blocks.
  [Default value = 50000]


  --extend_ab_by_adding_fragment N

  Extend sequences SA and SB: N = integer bigger than or equal to zero
  Extend sequences SA and SB by adding to them a fragment of length N 
  from the sides that are not supposed to be orthologous. It is equivalent 
  to add the non orthologous genes Co and Do to the sequences.
  Warning: This parameter is available only for table of synteny blocks.
  [Default value = 50000]


-------------------------------------------------------------------------------
File Formats:

  + Table orthologous genes:

  The table of orthologous genes must have just one 2 one orthologous genes 
  that can be found in the genomes GR and GO.

      - g1      = Name of the gene in the genome GR
      - c1      = Chromosome of the genome GR where g1 is located
      - inf1    = Start position of the gene g1 in the chromosome c1
      - sup1    = End position of the gene g1 in the chromosome c1
      - strand1 = Strand of the gene g1 (1 or -1)
      - g2      = Name of the gene in the genome GO
      - c2      = Chromosome of the genome GO where g2 is located
      - inf2    = Start position of the gene g2 in the chromosome c2
      - sup2    = End position of the gene g2 in the chromosome c2
      - strand2 = Strand of the gene g2 (1 or -1)

  Warning: The file must NOT have header with the column names.

  + Table orthologous synteny blocks:

  The table of orthologous synteny blocks must have just one 2 one 
  orthologous synteny blocks that can be found in the genomes GR and
  GO.

      - id      = Name of the synteny block
      - c1      = Chromosome of the genome GR where the block is located
      - inf1    = Start position of the block in the chromosome c1
      - sup1    = End position of the block in the chromosome c1
      - c2      = Chromosome of the genome GO where the block is located
      - inf2    = Start position of the block in the chromosome c2
      - sup2    = End position of the block in the chromosome c2
      - strand  = If 1, the synteny blocks of the genomes GR and GO are on
                  the same strand. Otherwise, they are on different strands.

  Warning: The file must NOT have header with the column names.

  + FASTA files:
  
  The fasta files must be named with the following format: 

      <chromosomename>.fasta

    Chromosome Name      File 
          1              1.fasta
          2A             2A.fasta
          X              X.fasta

  The file must have just one sequence that is related to the chromosome which
  is specified by the name of the file. If the FASTA file have more than one 
  sequence, the file will be ignored.

-------------------------------------------------------------------------------
SCRIPT OUTPUT:

  The script will write all results inside of the directory <outputdir>.
 
  Directories:

  + <outputdir>/alignments

    This directory will receive the results of the alignments of the sequences 
    SR vs SA and SR vs SB.

  + <outputdir>/dotplot

    This directory will receive the dotplot representation of the alignments of
    the sequences SR vs SA and SR vs SB.

  + <outputdir>/fasta

    This directory will receive the FASTA file of all sequences SR, SA and SB.

  + <outputdir>/segmentation

    This directory will receive all plots and results of the breakpoints that 
    were processed by the segmentation process.

  Files:

  + <outputdir>/NonRefinedBreakpoints.txt

    This file will receive all information about the breakpoints that were
    identified. It contains the data used by the segmentation process:

    id        - Breakpoint ID
    type      - Type of the breakpoint: inter or intra
    sRgeneA   - Name of the gene/block A in the sequence SR (genome GR)
    sRgeneB   - Name of the gene/block B in the sequence SR (genome GR)
    sRchr     - Chromosome of the genes/blocks A and B (genome GR)
    sRstrandA - Strand of the gene/block A (genome GR)
    sRstrandB - Strand of the gene/block B (genome GR)
    sRinf     - Inferior boundary of the sequence SR
    sRsup     - Superior boundary of the sequence SR
    sOgeneA   - Name of the gene/block A in the sequence SA (genome GO)
    sOgeneB   - Name of the gene/block B in the sequence SB (genome GO)
    sOchrA    - Chromosome of the gene/block A (genome GO)
    sOchrB    - Chromosome of the gene/block B (genome GO)
    sOstrandA - Strand of the gene/block A (genome GO)
    sOstrandB - Strand of the gene/block B (genome GO)
    sOinfA    - Inferior boundary of the sequence SA
    sOsupA    - Superior boundary of the sequence SA
    sOinfB    - Inferior boundary of the sequence SB
    sOsupB    - Superior boundary of the sequence SB
    bkpBegin  - Relative position of the breakpoint begin (related to sRinf)
    bkpEnd    - Relative position of the breakpoint end (related to sRinf)
    status    - Status of the breakpoint

    The status of the breakpoint can have one of the following values:
      1 : Valid breakpoint
     -2 : Sequence SR smaller than the allowed limit
     -3 : Sequence SA smaller than the allowed limit
     -4 : Sequence SB smaller than the allowed limit
     -5 : Sequence SA and SB smaller than the allowed limit
     -6 : Sequence SR bigger than the allowed limit

  + <outputdir>/segmentation.txt

    Final result of the segmentation process, this file will receive the 
    coordinates of the refined breakpoints:

    id        - Breakpoint id
    chr       - Chromosome where the breakpoint is located (Genome GR)
    oldbegin  - Old breakpoint begin position (in the chromosome sequence)
    oldend    - Old breakpoint end position (in the chromosome sequence)
    oldlength - Old breakpoint length
    newbegin  - New breakpoint begin position (in the chromosome sequence) 
                after the segmentation
    newend    - New breakpoint end position (in the chromosome sequence)
                after the segmentation
    newlength - New breakpoint length
    status    - Status of the breakpoint

    The status of the breakpoint can have one of the following values:
      1 : Segmentation result passed through the statistical test 
      0 : Segmentation result failed through the statistical test
     -1 : Segmentation was not performed because there are no hits on
          the alignments of the sequences SR vs SA and SR vs SB
     -2 : Segmentation was not performed because sequence SR is smaller
          than the allowed limit
     -3 : Segmentation was not performed because sequence SA is smaller 
          than the allowed limit
     -4 : Segmentation was not performed because sequence SB is smaller 
          than the allowed limit
     -5 : Segmentation was not performed because sequence SA and SB are 
          smaller than the allowed limit
     -6 : Segmentation was not performed because sequence SR is bigger 
          than the allowed limit
     -7 : Segmentation was not performed because R aborted the execution

End_print_usage

	if ( defined $error ) {
		print "\nERROR: $error\n\n";
	}

	exit(1);
}

###############################################################################
# This function prints an error message and aborts the execution
sub printErrorAndExit {
	my ($error) = @_;
	print STDERR "\nERROR: $error\n\n";
	exit(1);
}

###############################################################################
# This function validate the received parameters
sub validateParameters {

	my @parameters = @ARGV;
	my $n          = @parameters;

	if ( $n < 5 || ( $n > 5 && $n % 2 == 0 ) ) {
		printUsage("Invalid number of parameters.");
	}

	my $table   = fileExists( trim( $parameters[ $n - 5 ] ), \&printUsage );
	my $type    = trim( $parameters[ $n - 4 ] );
	my $dirGR   = directoryExists( trim( $parameters[ $n - 3 ] ), \&printUsage );
	my $dirGO   = directoryExists( trim( $parameters[ $n - 2 ] ), \&printUsage );
	my $workDir = directoryExists( trim( $parameters[ $n - 1 ] ), \&printUsage );

	if ( $type !~ /^[gGbB]$/ ) {
		printUsage("Parameter <type> : Invalid value");
	} else {
		$type =~ tr/gb/GB/;
	}

	# Mandatory parameters
	$PARAMETERS{"TABLE"}   = $table;
	$PARAMETERS{"TYPE"}    = $type;
	$PARAMETERS{"DIRGR"}   = $dirGR;
	$PARAMETERS{"DIRGO"}   = $dirGO;
	$PARAMETERS{"WORKDIR"} = $workDir;

	# Optional parameters : default values
	$PARAMETERS{"LASTZ_LEVEL"}    = $LASTZ_LEVEL;
	$PARAMETERS{"MAX_SR_LENGTH"}  = $MAXIMUM_SEQUENCE_R_SIZE;
	$PARAMETERS{"MAX_SAB_LENGTH"} = $MAXIMUM_SEQUENCE_AB_SIZE;
	$PARAMETERS{"EXTEND_BEFORE"}  = $EXTEND_BEFORE_VERIFY_LENGTH;

	if ( $type eq "G" ) {
		$PARAMETERS{"EXTEND"}     = $EXTEND_GENES;
		$PARAMETERS{"EXTEND_AB"}  = $EXTEND_AB_GENES;
		$PARAMETERS{"MIN_LENGTH"} = $MINIMUM_SEQUENCE_SIZE_GENE;
	} else {
		$PARAMETERS{"EXTEND"}     = $EXTEND_BLOCKS;
		$PARAMETERS{"EXTEND_AB"}  = $EXTEND_AB_BLOCKS;
		$PARAMETERS{"MIN_LENGTH"} = $MINIMUM_SEQUENCE_SIZE_BLOCKS;
	}

	for ( my $i = 0 ; $i < $n - 5 ; $i += 2 ) {
		my $parameter = trim( $parameters[$i] );
		my $value     = trim( $parameters[ $i + 1 ] );
		$parameter =~ tr/A-Z/a-z/;
		$value     =~ tr/a-z/A-Z/;

		if ( $parameter eq "--lastzlevel" ) {
			if ( $value !~ /^[123]$/ ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"LASTZ_LEVEL"} = $value;

		} elsif ( $parameter eq "--max_length_sr" ) {
			if ( $value !~ /^\d+$/ || $value == 0 ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"MAX_SR_LENGTH"} = $value;

		} elsif ( $parameter eq "--max_length_sab" ) {
			if ( $value !~ /^\d+$/ || $value == 0 ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"MAX_SAB_LENGTH"} = $value;

		} elsif ( $parameter eq "--minimum_length" ) {
			if ( $value !~ /^\d+$/ || $value == 0 ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"MIN_LENGTH"} = $value;

		} elsif ( $parameter eq "--extend_before" ) {
			if ( $value !~ /^[TF]$/ ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"EXTEND_BEFORE"} = $value;

		} elsif ( $parameter eq "--extend_by_adding_gene" ) {
			if ( $type eq "B" ) {
				printUsage(   "Parameter ${parameter} :\n"
							. "\tMust be set only when processing table of genes" );
			}
			if ( $value !~ /^[TF]$/ ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"EXTEND"} = $value;

		} elsif ( $parameter eq "--extend_ab_by_adding_gene" ) {
			if ( $type eq "B" ) {
				printUsage(   "Parameter ${parameter} :\n"
							. "\tMust be set only when processing table of genes" );
			}
			if ( $value !~ /^[TF]$/ ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"EXTEND_AB"} = $value;

		} elsif ( $parameter eq "--extend_by_adding_fragment" ) {
			if ( $type eq "G" ) {
				printUsage(   "Parameter ${parameter} :\n"
							. "\tMust be set only when processing table of blocks" );
			}
			if ( $value !~ /^\d+$/ ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"EXTEND"} = $value;

		} elsif ( $parameter eq "--extend_ab_by_adding_fragment" ) {
			if ( $type eq "G" ) {
				printUsage(   "Parameter ${parameter} :\n"
							. "\tMust be set only when processing table of blocks" );
			}
			if ( $value !~ /^\d+$/ ) {
				printUsage("Parameter ${parameter} : Invalid  value");
			}
			$PARAMETERS{"EXTEND_AB"} = $value;

		} else {
			printUsage("${parameter} : Unknown parameter");
		}
	}
}

###############################################################################
# This function verifies the environment (scripts, etc...)
sub verifyEnvironment {

	my $verifyDirScript   = fileExists( "${Bin}/core/verifyFASTAdirectory.pl",       \&printUsage );
	my $identifyScript    = fileExists( "${Bin}/core/identifyBreakpoints.pl",        \&printUsage );
	my $createFASTAScript = fileExists( "${Bin}/pipeline/createBreakpointFastas.pl", \&printUsage );
	my $alignScript       = fileExists( "${Bin}/pipeline/alignBreakpointRegions.pl", \&printUsage );
	my $dotplotScript = fileExists( "${Bin}/pipeline/dotplotBreakpointRegions.pl", \&printUsage );
	my $segmentationScript =
	  fileExists( "${Bin}/pipeline/breakpointSegmentation.pl", \&printUsage );

	$LASTZ_PATH = fileExists( $LASTZ_PATH, \&printUsage );

	if ( $LASTZ_LEVEL1_MATRIX !~ /^\s*default\s*$/ ) {
		$LASTZ_LEVEL1_MATRIX = fileExists( $LASTZ_LEVEL1_MATRIX, \&printUsage );
	}
	if ( $LASTZ_LEVEL2_MATRIX !~ /^\s*default\s*$/ ) {
		$LASTZ_LEVEL2_MATRIX = fileExists( $LASTZ_LEVEL2_MATRIX, \&printUsage );
	}
	if ( $LASTZ_LEVEL3_MATRIX !~ /^\s*default\s*$/ ) {
		$LASTZ_LEVEL3_MATRIX = fileExists( $LASTZ_LEVEL3_MATRIX, \&printUsage );
	}

	$ENVIRONMENT{"VERIFYDIR"}    = $verifyDirScript;
	$ENVIRONMENT{"IDENTIFY"}     = $identifyScript;
	$ENVIRONMENT{"CREATEFASTA"}  = $createFASTAScript;
	$ENVIRONMENT{"ALIGN"}        = $alignScript;
	$ENVIRONMENT{"DOTPLOT"}      = $dotplotScript;
	$ENVIRONMENT{"SEGMENTATION"} = $segmentationScript;

}

###############################################################################
# This function creates the directories where the script will write the results
sub createWorkDirectories {
	print "Creating work directories...\n";

	# Directory names
	my $workDir         = $PARAMETERS{"WORKDIR"};
	my $fastaDir        = "${workDir}/fasta";
	my $alignDir        = "${workDir}/alignments";
	my $dotplotDir      = "${workDir}/dotplot";
	my $segmentationDir = "${workDir}/segmentation";
	my $tmpDir          = "${workDir}/tmp";

	# Create the directories
	createDirectory($fastaDir);
	createDirectory($alignDir);
	createDirectory($dotplotDir);
	createDirectory($segmentationDir);
	createDirectory($tmpDir);

	# Verify if the directories were created
	if ( directoryExists($fastaDir) eq "" ) {
		printErrorAndExit("Could not create the directory ${fastaDir}");
	}
	if ( directoryExists($alignDir) eq "" ) {
		printErrorAndExit("Could not create the directory ${alignDir}");
	}
	if ( directoryExists($dotplotDir) eq "" ) {
		printErrorAndExit("Could not create the directory ${dotplotDir}");
	}
	if ( directoryExists($segmentationDir) eq "" ) {
		printErrorAndExit("Could not create the directory ${segmentationDir}");
	}
	if ( directoryExists($tmpDir) eq "" ) {
		printErrorAndExit("Could not create the directory ${tmpDir}");
	}

	# Clean the temporary directory
	deleteDirectoryContent("${tmpDir}");

	# Put the values in the parameters Hash
	$PARAMETERS{"FASTADIR"}        = $fastaDir;
	$PARAMETERS{"ALIGNDIR"}        = $alignDir;
	$PARAMETERS{"DOTPLOTDIR"}      = $dotplotDir;
	$PARAMETERS{"SEGMENTATIONDIR"} = $segmentationDir;
	$PARAMETERS{"TMPDIR"}          = $tmpDir;

}

###############################################################################
# This function validates the FASTA directories
sub validateFastaDirectories {

	my $dirGR  = $PARAMETERS{"DIRGR"};
	my $grFile = $PARAMETERS{"TMPDIR"} . "/" . rand() . time();

	my $dirGO  = $PARAMETERS{"DIRGO"};
	my $goFile = $PARAMETERS{"TMPDIR"} . "/" . rand() . time();

	print "Reading FASTA directory (Genome GR)...\n";
	system("perl $ENVIRONMENT{'VERIFYDIR'} ${dirGR} ${grFile}");

	print "Reading FASTA directory (Genome GO)...\n";
	system("perl $ENVIRONMENT{'VERIFYDIR'} ${dirGO} ${goFile}");

	if ( !-e $grFile || -s $grFile == 0 ) {
		printErrorAndExit("Could not find any FASTA file inside of ${dirGR}");
	}

	if ( !-e $goFile || -s $goFile == 0 ) {
		printErrorAndExit("Could not find any FASTA file inside of ${dirGO}");
	}

	$PARAMETERS{"GRCHR"} = $grFile;
	$PARAMETERS{"GOCHR"} = $goFile;
}

###############################################################################
# This function validates the input table
sub validateInputTable {

	my $table    = $PARAMETERS{"TABLE"};
	my $type     = $PARAMETERS{"TYPE"};
	my $validate = $PARAMETERS{"TMPDIR"} . "/" . rand() . time();

	print "Verifying the input table\n";

	# Read the chromosome list of the genome GR
	my %chrGR = readChromosomeFile( $PARAMETERS{"GRCHR"} );

	# Read the chromosome list of the genome GO
	my %chrGO = readChromosomeFile( $PARAMETERS{"GOCHR"} );

	if ( $type eq "G" ) {

		my %genes1;
		my %genes2;

		# Table of orthologous genes
		my $nLine = 0;
		open( IN, "${table}" ) or die printErrorAndExit("Could not read ${table}");
		open( OUT, ">${validate}" )
		  or die printErrorAndExit("Could not write temporary file ${validate}");
		while ( my $line = <IN> ) {
			$nLine++;
			$line = trim($line);
			my @columns = split( "\t", $line );
			my $nColumns = @columns;
			if ( $nColumns == 10 ) {
				my $g1      = trim( $columns[0] );
				my $chr1    = trim( $columns[1] );
				my $inf1    = trim( $columns[2] );
				my $sup1    = trim( $columns[3] );
				my $strand1 = trim( $columns[4] );
				my $g2      = trim( $columns[5] );
				my $chr2    = trim( $columns[6] );
				my $inf2    = trim( $columns[7] );
				my $sup2    = trim( $columns[8] );
				my $strand2 = trim( $columns[9] );

				$chr1 =~ tr/a-z/A-Z/;
				$chr2 =~ tr/a-z/A-Z/;

				if ( $inf1 > $sup1 ) {
					my $aux = $inf1;
					$inf1 = $sup1;
					$sup1 = $aux;
				}

				if ( $inf2 > $sup2 ) {
					my $aux = $inf2;
					$inf2 = $sup2;
					$sup2 = $aux;
				}

				my $ok = 1;
				if ( $strand1 !~ /^[\+\-]1*$/ && $strand1 !~ /^[\+\-]?1$/ ) {
					print "Warning: File ${table} - Line ${nLine} - Column 5:\n";
					print "\tInvalid strand value at column 5. Line was ignored.\n";
					$ok = 0;
				}

				if ( $strand1 =~ /^\+$/ ) {
					$strand1 = 1;
				} elsif ( $strand1 =~ /^\-$/ ) {
					$strand1 = -1;
				}

				if ( $ok == 1 && $strand2 !~ /^[\+\-]1*$/ && $strand2 !~ /^[\+\-]?1$/ ) {
					print "Warning: File ${table} - Line ${nLine} - Column 10:\n";
					print "\tInvalid strand value at column 10. Line was ignored.\n";
					$ok = 0;
				}

				if ( $strand2 =~ /^\+$/ ) {
					$strand2 = 1;
				} elsif ( $strand2 =~ /^\-$/ ) {
					$strand2 = -1;
				}

				if ( $ok == 1 && !defined $chrGR{$chr1} ) {
					print "Warning: File ${table} - Line ${nLine} - Column 2:\n";
					print "\tChromosome ${chr1} does not have FASTA. Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 && !defined $chrGO{$chr2} ) {
					print "Warning: File ${table} - Line ${nLine} - Column 7:\n";
					print "\tChromosome ${chr2} does not have FASTA. Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 && defined $genes1{$g1} ) {
					print "Warning: File ${table} - Line ${nLine} - Column 1:\n";
					print "\tGene ${g1} appears more than one time. Line was ignored\n";
					$ok = 0;
				}
				$genes1{$g1} = 1;

				if ( $ok == 1 && defined $genes2{$g2} ) {
					print "Warning: File ${table} - Line ${nLine} - Column 6:\n";
					print "\tGene ${g2} appears more than one time. Line was ignored\n";
					$ok = 0;
				}
				$genes2{$g2} = 1;

				if ( $ok == 1 && ( $inf1 < 1 || $sup1 > $chrGR{$chr1} ) ) {
					print "Warning: File ${table} - Line ${nLine} - Columns 3 and 4:\n";
					print "\tInvalid coordinates for the gene ${g1}. Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 && ( $inf2 < 1 || $sup2 > $chrGO{$chr2} ) ) {
					print "Warning: File ${table} - Line ${nLine} - Columns 8 and 9:\n";
					print "\tInvalid coordinates for the gene ${g2}. Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 ) {
					my @a =
					  ( $g1, $chr1, $inf1, $sup1, $strand1, $g2, $chr2, $inf2, $sup2, $strand2 );
					print OUT join( "\t", @a ) . "\n";
				}

			} else {
				print "Warning: File ${table} - Line ${nLine}:\n";
				print "\tInvalid number of columns. Line was ignored.\n";
			}

		}
		close(OUT);
		close(IN);

		$PARAMETERS{"VALIDATEDTABLE"} = $validate;

	} else {

		my %blocks;

		# Synteny blocks
		my $nLine = 0;
		open( IN, "${table}" ) or die printErrorAndExit("Could not read ${table}");
		open( OUT, ">${validate}" )
		  or die printErrorAndExit("Could not write temporary file ${validate}");
		while ( my $line = <IN> ) {
			$nLine++;
			$line = trim($line);
			my @columns = split( "\t", $line );
			my $nColumns = @columns;
			if ( $nColumns == 8 ) {
				my $blockId = trim( $columns[0] );
				my $chr1    = trim( $columns[1] );
				my $inf1    = trim( $columns[2] );
				my $sup1    = trim( $columns[3] );
				my $chr2    = trim( $columns[4] );
				my $inf2    = trim( $columns[5] );
				my $sup2    = trim( $columns[6] );
				my $strand  = trim( $columns[7] );

				$chr1 =~ tr/a-z/A-Z/;
				$chr2 =~ tr/a-z/A-Z/;

				if ( $inf1 > $sup1 ) {
					my $aux = $inf1;
					$inf1 = $sup1;
					$sup1 = $aux;
				}

				if ( $inf2 > $sup2 ) {
					my $aux = $inf2;
					$inf2 = $sup2;
					$sup2 = $aux;
				}

				my $ok = 1;
				if ( $strand !~ /^[\+\-]1*$/ && $strand !~ /^[\+\-]?1$/ ) {
					print "Warning: File ${table} - Line ${nLine} - Column 5:\n";
					print "\tInvalid strand value at column 8. Line was ignored.\n";
					$ok = 0;
				}

				if ( $strand =~ /^\+$/ ) {
					$strand = 1;
				} elsif ( $strand =~ /^\-$/ ) {
					$strand = -1;
				}

				if ( $ok == 1 && !defined $chrGR{$chr1} ) {
					print "Warning: File ${table} - Line ${nLine} - Column 2:\n";
					print "\tChromosome ${chr1} does not have FASTA. Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 && !defined $chrGO{$chr2} ) {
					print "Warning: File ${table} - Line ${nLine} - Column 5:\n";
					print "\tChromosome ${chr2} does not have FASTA. Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 && defined $blocks{$blockId} ) {
					print "Warning: File ${table} - Line ${nLine} - Column 1:\n";
					print "\tBlock ${blockId} appears more than one time. Line was ignored\n";
					$ok = 0;
				}
				$blocks{$blockId} = 1;

				if ( $ok == 1 && ( $inf1 < 1 || $sup1 > $chrGR{$chr1} ) ) {
					print "Warning: File ${table} - Line ${nLine} - Columns 3 and 4:\n";
					print "\tInvalid coordinates for the block ${blockId}. ";
					print "Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 && ( $inf2 < 1 || $sup2 > $chrGO{$chr2} ) ) {
					print "Warning: File ${table} - Line ${nLine} - Columns 6 and 7:\n";
					print "\tInvalid coordinates for the block ${blockId}. ";
					print "Line was ignored\n";
					$ok = 0;
				}

				if ( $ok == 1 ) {
					my $strand1 = 1;
					my $strand2 = 1;
					if ( $strand =~ /^\-/ ) {
						$strand2 = -1;
					}
					my @a = (
							  "${blockId}.R", $chr1, $inf1, $sup1, $strand1, "${blockId}.O", $chr2,
							  $inf2, $sup2, $strand2
					);
					print OUT join( "\t", @a ) . "\n";
				}

			} else {
				print "Warning: File ${table} - Line ${nLine}:\n";
				print "\tInvalid number of columns. Line was ignored.\n";
			}

		}
		close(OUT);
		close(IN);

		$PARAMETERS{"VALIDATEDTABLE"} = $validate;
	}
}

###############################################################################
# This function reads the chromosome file to find the end position of each
# chromosome
sub readChromosomeFile {
	my ($file) = @_;
	my %toReturn;
	open( IN, "${file}" )
	  or die printErrorAndExit("Could not read temporary file ${file}");
	while ( my $line = <IN> ) {
		chomp($line);
		my ( $chr, $begin, $end, $file ) = split( "\t", $line );
		$toReturn{$chr} = $end;
	}
	close(IN);
	return %toReturn;
}

###############################################################################
# This function calls the script that identify the breakpoint regions
sub identifyBreakpoints {

	my $script       = $ENVIRONMENT{"IDENTIFY"};
	my $table        = $PARAMETERS{"VALIDATEDTABLE"};
	my $type         = $PARAMETERS{"TYPE"};
	my $grFile       = $PARAMETERS{"GRCHR"};
	my $goFile       = $PARAMETERS{"GOCHR"};
	my $workDir      = $PARAMETERS{"WORKDIR"};
	my $tmpDir       = $PARAMETERS{"TMPDIR"};
	my $breakpoints  = $workDir . "/NonRefinedBreakpoints.txt";
	my $extend       = $PARAMETERS{"EXTEND"};
	my $extendAB     = $PARAMETERS{"EXTEND_AB"};
	my $minLength    = $PARAMETERS{"MIN_LENGTH"};
	my $maxSR        = $PARAMETERS{"MAX_SR_LENGTH"};
	my $maxSAB       = $PARAMETERS{"MAX_SAB_LENGTH"};
	my $extendBefore = $PARAMETERS{"EXTEND_BEFORE"};

	print "Identyfing breakpoints...\n";
	if (-s $table == 0) {
		printErrorAndExit("Breakpoint identification step failed. Input table is empty.");
	}
	my $cmd = "perl ${script} -i ${table} -t ${type} -o ${breakpoints} ";
	$cmd .= "-R ${grFile} -O ${goFile} -T ${tmpDir} -M ${minLength} ";
	$cmd .= "-l ${maxSR} -L ${maxSAB} -v ${extendBefore} ";
	if ( $type eq "G" ) {
		$cmd .= "-e ${extend} -E ${extendAB} ";
	} else {
		$cmd .= "-b ${extend} -B ${extendAB} ";
	}

	system($cmd);
	deleteFile("${breakpoints}.log.Rout");

	$breakpoints = fileExists($breakpoints);
	if ( $breakpoints eq "" ) {
		printErrorAndExit("Breakpoint identification step failed.");
	}
	$PARAMETERS{"BREAKPOINTTABLE"} = $breakpoints;

}

###############################################################################
# This function calls the script that create the breakpoint regions
sub createBreakpointFASTAs {

	my $breakpoints = $PARAMETERS{"BREAKPOINTTABLE"};
	my $fastaDir    = $PARAMETERS{"FASTADIR"};
	my $grFile      = $PARAMETERS{"GRCHR"};
	my $goFile      = $PARAMETERS{"GOCHR"};
	my $script      = $ENVIRONMENT{"CREATEFASTA"};

	print "Creating sequences SR, SA and SB (FASTA files)...\n";
	my $cmd = "perl ${script} ${breakpoints} ${grFile} ${goFile} ${fastaDir}";
	system($cmd);

}

###############################################################################
# This function calls the script which aligns the breakpoint sequences
sub alignBreakpoints {

	my $script     = $ENVIRONMENT{"ALIGN"};
	my $fastaDir   = $PARAMETERS{"FASTADIR"};
	my $alignDir   = $PARAMETERS{"ALIGNDIR"};
	my $lastzlevel = $PARAMETERS{"LASTZ_LEVEL"};
	my $tmpDir     = $PARAMETERS{"TMPDIR"};

	my $lastzParameters = $LASTZ_LEVEL2_PARAMETERS;
	my $lastzMatrix     = $LASTZ_LEVEL2_MATRIX;
	if ( $lastzlevel == 1 ) {
		$lastzParameters = $LASTZ_LEVEL1_PARAMETERS;
		$lastzMatrix     = $LASTZ_LEVEL1_MATRIX;
	} elsif ( $lastzlevel == 3 ) {
		$lastzParameters = $LASTZ_LEVEL3_PARAMETERS;
		$lastzMatrix     = $LASTZ_LEVEL3_MATRIX;
	}

	print "Aligning breakpoint sequences (SR vs SA and SR vs SB)...\n";
	my $cmd = "perl ${script} ${fastaDir} ${alignDir} ";
	$cmd .= "\"${lastzParameters}\" ${lastzMatrix} ${LASTZ_PATH} ${tmpDir}";
	system($cmd);

}

###############################################################################
# This function calls the script which create dotplot for the breakpoints
sub dotplotBreakpoints {

	my $breakpoints = $PARAMETERS{"BREAKPOINTTABLE"};
	my $script      = $ENVIRONMENT{"DOTPLOT"};
	my $fastaDir    = $PARAMETERS{"FASTADIR"};
	my $alignDir    = $PARAMETERS{"ALIGNDIR"};
	my $dotplotDir  = $PARAMETERS{"DOTPLOTDIR"};
	my $tmpDir      = $PARAMETERS{"TMPDIR"};

	print "Creating dotplots...\n";
	my $cmd = "perl ${script} ${breakpoints} ${fastaDir} ${alignDir} ";
	$cmd .= "${dotplotDir} ${tmpDir}";
	system($cmd);

}

###############################################################################
# This function calls the script which refines the breakpoints
sub refineBreakpoints {

	my $breakpoints     = $PARAMETERS{"BREAKPOINTTABLE"};
	my $script          = $ENVIRONMENT{"SEGMENTATION"};
	my $fastaDir        = $PARAMETERS{"FASTADIR"};
	my $alignDir        = $PARAMETERS{"ALIGNDIR"};
	my $segmentationDir = $PARAMETERS{"SEGMENTATIONDIR"};
	my $tmpDir          = $PARAMETERS{"TMPDIR"};
	my $workDir         = $PARAMETERS{"WORKDIR"};

	# Run the script
	print "Refining breakpoint regions (segmentation)...\n";
	my $cmd = "perl ${script} ${breakpoints} ${fastaDir} ${alignDir} ";
	$cmd .= "${segmentationDir} ${tmpDir}";
	system($cmd);

	# Create a single file with the result of the segmentation of all breakpoints
	my $resultFile = "${workDir}/segmentation.txt";
	deleteFile($resultFile);

	my @files = `ls -1 ${segmentationDir}`;
	chomp(@files);
	@files = sort @files;
	foreach my $file (@files) {
		if ( $file =~ /^Breakpoint_\d+\.txt$/ ) {
			`cat ${segmentationDir}/${file} >> ${resultFile}`;
		}
	}

	# Clean the temporary directory
	deleteDirectoryContent("${tmpDir}");

}
