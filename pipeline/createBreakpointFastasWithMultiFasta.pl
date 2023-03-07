#!/usr/bin/perl -w

###############################################################################
# Copyright (C) 2023 Inria                        #
# 									      #
# Contributors : Christian Baudet, Claire Lemaitre			      #
# 									      #
# Contact : christian.baudet@gmail.com, claire.lemaitre@inria.fr	      #
# 									      #
# This file is part of Cassis.						      #
# 									      #
# Cassis is free software: you can redistribute it and/or modify	      #
# it under the terms of the GNU General Public License as published by	      #
#  the Free Software Foundation, either version 3 of the License, or	      #
# (at your option) any later version.					      #
# 									      #
# Cassis is distributed in the hope that it will be useful,		      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of	      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		      #
# GNU General Public License for more details.				      #
# 									      #
# You should have received a copy of the GNU General Public License	      #
# along with Cassis.  If not, see <http://www.gnu.org/licenses/>.             #
###############################################################################

use strict;
use File::Spec;
use FindBin qw($Bin);
require "${Bin}/../core/util.pm";

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################

# Verify the environment
my %ENVIRONMENT;
verifyEnvironment();

# Receive and process the parameters
my %PARAMETERS;
validateParameters();

# Create the SR fasta files
createSRfastas();

# Create the SA and SB fasta files
createSABfastas();

exit(0);

###############################################################################
## AUXILIARY SUB-ROUTINES #####################################################
###############################################################################

###############################################################################
# This function prints the script usage mode
sub printUsage {
  
  my ($message) = @_;
  
  if (defined $message) {
    print STDERR "\nERROR: $message\n";
  }
  
  print << "End_Print_Usage";

Usage:
  perl createBreakpointFastasWithMultiFasta.pl <table> <fastaGR> <fastaGO> <outputdir>

-------------------------------------------------------------------------------
Mandatory parameters:

  <table> Path for the file that contains the list of breakpoints.

  <fastaGR> Path for the fasta file that contains the chromosome sequences of the whole
          genome GR.

  <fastaGO> Path for the fasta file that contains the chromosome sequences of the whole
          genome GO.

  <outputdir> Path for the directory where the script will write the FASTA
              files of the breakpoint regions.

-------------------------------------------------------------------------------
File formats:

  + Breakpoint table file:

  The breakpoint table file has the information about the breakpoints that were
  verified between two genomes GR and GO. It has the following columns:

  + id        - Breakpoint ID
  + type      - Type of the breakpoint: inter or intra
  + sRgeneA   - Name of the gene/block A in the sequence SR (genome GR)
  + sRgeneB   - Name of the gene/block B in the sequence SR (genome GR)
  + sRchr     - Chromosome of the genes/blocks A and B (genome GR)
  + sRstrandA - Strand of the gene/block A (genome GR)
  + sRstrandB - Strand of the gene/block B (genome GR)
  + sRinf     - Inferior boundary of the sequence SR
  + sRsup     - Superior boundary of the sequence SR
  + sOgeneA   - Name of the gene/block A in the sequence SA (genome GO)
  + sOgeneB   - Name of the gene/block B in the sequence SB (genome GO)
  + sOchrA    - Chromosome of the gene/block A (genome GO)
  + sOchrB    - Chromosome of the gene/block B (genome GO)
  + sOstrandA - Strand of the gene/block A (genome GO)
  + sOstrandB - Strand of the gene/block B (genome GO)
  + sOinfA    - Inferior boundary of the sequence SA
  + sOsupA    - Superior boundary of the sequence SA
  + sOinfB    - Inferior boundary of the sequence SB
  + sOsupB    - Superior boundary of the sequence SB
  + bkpBegin  - Relative position of the breakpoint begin (related to sRinf)
  + bkpEnd    - Relative position of the breakpoint end (related to sRinf)
  + status    - Status of the breakpoint

  The status of the breakpoint can have one of the following values:
     1 : Valid breakpoint
    -2 : Sequence SR smaller than the allowed limit
    -3 : Sequence SA smaller than the allowed limit
    -4 : Sequence SB smaller than the allowed limit
    -5 : Sequence SA and SB smaller than the allowed limit
    -6 : Sequence SR bigger than the allowed limit

  Warning: The file must NOT have header with the column names.

  + multi-fasta files <fastaGR> <fastaGO>

-------------------------------------------------------------------------------
Script output:

The script will read the list of breakpoints of the file <table> and produce 
all sequences SR, SA and SB. It will write inside of the directory <outputdir> 
the sequences (SR, SA and SB) which are present in the chromosomes that are 
listed in the files <chrGR> and <chrGO>.

The script will create FASTA files, inside of the directory <outputdir>, that 
have the following name pattern: Breakpoint_[N]_S[X].fasta, where [N] is the 
number of the breakpoint and [X] is one of the letters R, A or, B.

End_Print_Usage
  
  if (defined $message) {
    print "\nERROR: $message\n\n";
  }
  
  exit(1);
}

###############################################################################
# This function validate the received parameters
sub validateParameters {

  my @parameters = @ARGV;
  my $nParameters = @parameters;
  
  if ($nParameters != 4) {
    printUsage("Invalid number of parameters.");
  }
  
  my $table     = fileExists(trim($parameters[0]), \&printUsage);
  my $fastaGR     = fileExists(trim($parameters[1]), \&printUsage);
  my $fastaGO     = fileExists(trim($parameters[2]), \&printUsage);
  my $outputdir = directoryExists(trim($parameters[3]), \&printUsage);

  $PARAMETERS{"TABLE"}     = $table;
  $PARAMETERS{"FASTAGR"}     = $fastaGR;
  $PARAMETERS{"FASTAGO"}     = $fastaGO;
  $PARAMETERS{"OUTPUTDIR"} = $outputdir;

}

###############################################################################
# This function verifies with the core script is ok
sub verifyEnvironment {
  my $corescript = "${Bin}/../core/extractBreakpointRegionsWithMultiFasta.pl";
  if (! -e $corescript) {
    printUsage("Missing script ${corescript}");
  }
  $ENVIRONMENT{"CORESCRIPT"} = $corescript;
}

###############################################################################
# This function creates the fasta files of the sequences SR
sub createSRfastas {

  my $table     = $PARAMETERS{"TABLE"};
  my $fastaGR     = $PARAMETERS{"FASTAGR"};
  my $outputdir = $PARAMETERS{"OUTPUTDIR"};
  my $script    = $ENVIRONMENT{"CORESCRIPT"};

  my $cmd = "perl ${script} ${table} ${fastaGR} ${outputdir} 1";
  print "Creating sequences SR from Genome GR\n";
  system("$cmd");

}

###############################################################################
# This function creates the fasta files of the sequences SA and SB
sub createSABfastas {

  my $table     = $PARAMETERS{"TABLE"};
  my $fastaGO     = $PARAMETERS{"FASTAGO"};
  my $outputdir = $PARAMETERS{"OUTPUTDIR"};
  my $script    = $ENVIRONMENT{"CORESCRIPT"};

    my $cmd = "perl ${script} ${table} ${fastaGO} ${outputdir} 2";
    print "Creating sequences SA and SB from Genome GO\n";
    system("$cmd");

}
