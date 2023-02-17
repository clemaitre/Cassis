#!/usr/bin/perl -w

###############################################################################
# Copyright (C) 2010 Université Claude Bernard Lyon 1                         #
# 									      #
# Contributors : Christian Baudet, Claire Lemaitre			      #
# 									      #
# Contact : christian.baudet@gmail.com, claire.lemaitre@gmail.com	      #
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
require "${Bin}/util.pm";

###############################################################################
# Constant: R command line
my $R_CMD = "R CMD BATCH --no-restore --no-save";

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################

# Verify the environment
my %ENVIRONMENT;
verifyEnvironment();

# Receive and process the parameters
my %PARAMETERS;
validateParameters();

# Generate the dotplot
dotPlotBreakpoint();

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
  perl dotplotBreakpoint.pl <number> <table> <lastzSA> <lastzSB> <fastaSR> <pngFile> [tmpdir]

-------------------------------------------------------------------------------
Mandatory parameters:

  <number>  Number of the breakpoint.

  <table>   Path for the file that contains the list of breakpoints.

  <lastzSA> Path for the file that has the results of the alignment of SR
            against SA.

  <lastzSB> Path for the file that has the results of the alignment of SR
            against SB.

  <fastaSR> Path for the FASTA file of the sequence SR.

  <pngfile> Path for the file that will receive the dotplot.

-------------------------------------------------------------------------------
Optional parameter:

  [tmpdir] Path for the temporary directory. /tmp/ is the default temporary 
           directory.

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

  + LASTZ alignment result file:

  The LASTZ alignment result file has the information about the hits of the
  aligment of the sequence SR against SA (or SB). It has the following fields:

  + score - Hit score
  + b1    - Begin position of the hit in the first sequence (SR)
  + b2    - Begin position of the hit in the second sequence (SA or SB)
  + e1    - End position of the hit in the first sequence (SR)
  + e2    - End position of the hit in the second sequence (SA or SB)

  Warning: The file must NOT have header with the column names.

-------------------------------------------------------------------------------
Script output:

The script create a dotplot that show the distribution of the hits of the
alignment between the sequences SR vs SA and SR vs SB. 

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
  
  if ($nParameters < 6 || $nParameters > 7) {
    printUsage("Invalid number of parameters.");
  }
  
  my $number   = trim($parameters[0]);
  my $table    = fileExists(trim($parameters[1]), \&printUsage);
  my $lastzSA  = fileExists(trim($parameters[2]), \&printUsage);
  my $lastzSB  = fileExists(trim($parameters[3]), \&printUsage);
  my $fastaSR  = fileExists(trim($parameters[4]), \&printUsage);
  my $pngFile  = trim($parameters[5]);
  my $tmpDir   = "/tmp";
  if ($nParameters == 7) {
    $tmpDir = directoryExists(trim($parameters[6]), \&printUsage);
  }
  
  if ($number !~ /^(\d+)$/) {
    printUsage("Invalid value for the parameter <number>");
  }

  $PARAMETERS{"NUMBER"}   = $number;
  $PARAMETERS{"TABLE"}    = $table;
  $PARAMETERS{"LASTZSA"}  = $lastzSA;
  $PARAMETERS{"LASTZSB"}  = $lastzSB;
  $PARAMETERS{"FASTASR"}  = $fastaSR;
  $PARAMETERS{"PNG"}      = $pngFile;
  $PARAMETERS{"TMPDIR"}   = $tmpDir;

}

###############################################################################
# This function verifies with the R scripts are ok
sub verifyEnvironment {

  my $dotplotScript        = fileExists("${Bin}/R/dotplot.R", \&printUsage);
  my $utilScript           = fileExists("${Bin}/R/util.R", \&printUsage);
  my $intervalsScript      = fileExists("${Bin}/R/intervals.R", \&printUsage);
  my $maskedIntervalScript = fileExists("${Bin}/getMaskedIntervals.pl", \&printUsage);
  $ENVIRONMENT{"DOTPLOT"}  = $dotplotScript;
  $ENVIRONMENT{"UTIL"}     = $utilScript;
  $ENVIRONMENT{"INTERVAL"} = $intervalsScript;
  $ENVIRONMENT{"MASKED"}   = $maskedIntervalScript;

}

###############################################################################
# This function generates the dotplot for the breakpoint region
sub dotPlotBreakpoint {

  my $number       = $PARAMETERS{"NUMBER"};
  my $table        = $PARAMETERS{"TABLE"};
  my $lastzSA      = $PARAMETERS{"LASTZSA"};
  my $lastzSB      = $PARAMETERS{"LASTZSB"};
  my $fastaSR      = $PARAMETERS{"FASTASR"};
  my $pngFile      = $PARAMETERS{"PNG"};
  my $tmpDir       = $PARAMETERS{"TMPDIR"};
  my $maskedScript = $ENVIRONMENT{"MASKED"};

  my $maskedIntervals = "${tmpDir}/".rand().time();
  my $cmd = "perl ${maskedScript} ${fastaSR} ${maskedIntervals}";
  system($cmd);

  if (-s $maskedIntervals == 0) {
    $maskedIntervals = "";
  }
  if (-s $lastzSA == 0) {
    $lastzSA = "";
  }
  if (-s $lastzSB == 0) {
    $lastzSB = "";
  }

  runScriptR($number, $lastzSA, $lastzSB, $table, $pngFile, $maskedIntervals);
  
  deleteFile($maskedIntervals);

}

###############################################################################
# This function runs the script R with the received parameters
sub runScriptR {

  my ($number, $lastzSA, $lastzSB, $table, $pngFile, $maskedIntervals) = @_;

  my $tmpDir          = $PARAMETERS{"TMPDIR"};
  my $dotplotScript   = $ENVIRONMENT{"DOTPLOT"};
  my $utilScript      = $ENVIRONMENT{"UTIL"};
  my $intervalsScript = $ENVIRONMENT{"INTERVAL"};
  my $scriptFile      = "${tmpDir}/".rand().time();
  my $png             = "${tmpDir}/".rand().time();

  my $n = int $number;

  # Write the script file
  open(OUT, ">${scriptFile}");

  print OUT "source(\"${dotplotScript}\")\n";
  print OUT "source(\"${utilScript}\")\n";
  print OUT "source(\"${intervalsScript}\")\n";
  print OUT "breakpointTable = readBreakpointsFile(\"${table}\", hasHeader=F)\n";

  if ($lastzSA ne "") {
    print OUT "alignmentSRSA = readLastzOutputFile(\"${lastzSA}\", hasHeader=F)\n";
  } else {
    print OUT "alignmentSRSA = data.frame()\n";
  }

  if ($lastzSB ne "") {
    print OUT "alignmentSRSB = readLastzOutputFile(\"${lastzSB}\", hasHeader=F)\n";
  } else {
    print OUT "alignmentSRSB = data.frame()\n";
  }

  if ($maskedIntervals ne "") {
    print OUT "maskedIntervals = readMaskedIntervalsFile(\"${maskedIntervals}\", hasHeader=F)\n";
  } else {
    print OUT "maskedIntervals = NULL\n";
  }

  print OUT "bitmap(\"${png}\", type=\"png16m\", height=8, width=7, res=100)\n";

  print OUT << "EndScript";

dotplotBreakpointRegion(breakpointId    = ${n},
                        alignmentSRSA   = alignmentSRSA,
                        alignmentSRSB   = alignmentSRSB,
                        breakpointTable = breakpointTable,
                        maskedIntervals = maskedIntervals,
                        geneList        = NULL)
dev.off()

EndScript

  close(OUT);

  # Execute R with the script
  my $logR = "${scriptFile}.log.Rout";
  system("$R_CMD ${scriptFile} ${logR}");
  if (! -e $png) {
    print STDERR "Warning: Breakpoint region number ${number} - Dotplot failed\n";
  } else {
    `mv ${png} ${pngFile}`;
  }
  deleteFile($scriptFile);
  deleteFile($logR);
  deleteFile($png);
}
