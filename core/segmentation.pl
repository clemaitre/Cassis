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

# Perform the segmentation
segmentation();

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
  perl segmentation.pl <number> <table> <lastzSA> <lastzSB> <fastaSR> <txtFile> <pngFile> [tmpdir]

-------------------------------------------------------------------------------
Mandatory parameters:

  <number>  Number of the breakpoint.

  <table>   Path for the file that contains the list of breakpoints.

  <lastzSA> Path for the file that has the results of the alignment of SR
            against SA.

  <lastzSB> Path for the file that has the results of the alignment of SR
            against SB.

  <fastaSR> Path for the FASTA file of the sequence SR.

  <txtfile> Path for the file that will receive the segmentation result.

  <pngfile> Path for the file that will receive the segmentation plot.

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

  The script will create a plot (<pngFile>) that is a graphical representation
  of the segmentation. It will generate also a TSV file (<txtFile>) that have
  the following columns:

  + id        - Breakpoint id
  + chr       - Chromosome where the breakpoint is located (Genome GR)
  + oldbegin  - Old breakpoint begin position (in the chromosome sequence)
  + oldend    - Old breakpoint end position (in the chromosome sequence)
  + oldlength - Old breakpoint length
  + newbegin  - New breakpoint begin position (in the chromosome sequence) 
  +             after the segmentation
  + newend    - New breakpoint end position (in the chromosome sequence)
  +             after the segmentation
  + newlength - New breakpoint length
  + status    - Status of the breakpoint

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
  
  if ($nParameters < 7 || $nParameters > 8) {
    printUsage("Invalid number of parameters.");
  }
  
  my $number   = trim($parameters[0]);
  my $table    = fileExists(trim($parameters[1]), \&printUsage);
  my $lastzSA  = fileExists(trim($parameters[2]), \&printUsage);
  my $lastzSB  = fileExists(trim($parameters[3]), \&printUsage);
  my $fastaSR  = fileExists(trim($parameters[4]), \&printUsage);
  my $txtFile  = trim($parameters[5]);
  my $pngFile  = trim($parameters[6]);
  my $tmpDir   = "/tmp";
  if ($nParameters == 8) {
    $tmpDir = directoryExists(trim($parameters[7]), \&printUsage);
  }
  
  if ($number !~ /^(\d+)$/) {
    printUsage("Invalid value for the parameter <number>");
  }

  $PARAMETERS{"NUMBER"}   = $number;
  $PARAMETERS{"TABLE"}    = $table;
  $PARAMETERS{"LASTZSA"}  = $lastzSA;
  $PARAMETERS{"LASTZSB"}  = $lastzSB;
  $PARAMETERS{"FASTASR"}  = $fastaSR;
  $PARAMETERS{"TXT"}      = $txtFile;
  $PARAMETERS{"PNG"}      = $pngFile;
  $PARAMETERS{"TMPDIR"}   = $tmpDir;

}

###############################################################################
# This function verifies with the R scripts are ok
sub verifyEnvironment {

  my $segmentationScript   = fileExists("${Bin}/R/segmentation.R", \&printUsage);
  my $utilScript           = fileExists("${Bin}/R/util.R", \&printUsage);
  my $intervalsScript      = fileExists("${Bin}/R/intervals.R", \&printUsage);
  my $maskedIntervalScript = fileExists("${Bin}/getMaskedIntervals.pl", \&printUsage);
  $ENVIRONMENT{"SEGMENTATION"} = $segmentationScript;
  $ENVIRONMENT{"UTIL"}         = $utilScript;
  $ENVIRONMENT{"INTERVAL"}     = $intervalsScript;
  $ENVIRONMENT{"MASKED"}       = $maskedIntervalScript;

}

###############################################################################
# This function calculates the segmentation of the breakpoint region
sub segmentation {
  
  my $number       = $PARAMETERS{"NUMBER"};
  my $table        = $PARAMETERS{"TABLE"};
  my $lastzSA      = $PARAMETERS{"LASTZSA"};
  my $lastzSB      = $PARAMETERS{"LASTZSB"};
  my $fastaSR      = $PARAMETERS{"FASTASR"};
  my $txtFile      = $PARAMETERS{"TXT"};
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

  runScriptR($number, $lastzSA, $lastzSB, $table, $pngFile,
	     $txtFile, $maskedIntervals);
  
  deleteFile($maskedIntervals);

}

###############################################################################
# This function runs the script R with the received parameters
sub runScriptR {

  my ($number, $lastzSA, $lastzSB, $table, $pngFile,
      $txtFile, $maskedIntervals) = @_;

  my $tmpDir             = $PARAMETERS{"TMPDIR"};
  my $segmentationScript = $ENVIRONMENT{"SEGMENTATION"};
  my $utilScript         = $ENVIRONMENT{"UTIL"};
  my $intervalsScript    = $ENVIRONMENT{"INTERVAL"};
  my $scriptFile         = "${tmpDir}/".rand().time();
  my $png                = "${tmpDir}/".rand().time();
  my $seg                = "${tmpDir}/".rand().time();

  my $n = int $number;

  # Write the script file
  open(OUT, ">${scriptFile}");

  print OUT "source(\"${segmentationScript}\")\n";
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

data = segmentAndPlotABreak(breakpointId    = ${n},
                            alignmentSRSA   = alignmentSRSA,
                            alignmentSRSB   = alignmentSRSB,
                            breakpointTable = breakpointTable,
                            maskedIntervals = maskedIntervals,
                            toplot          = T)

output = data.frame(id           = data\$id,
                    chrR         = data\$chrR,
                    oldBkpBegin  = data\$oldBkpBegin,
                    oldBkpEnd    = data\$oldBkpEnd,
                    oldBkpLength = data\$oldBkpLength,
                    bkpBegin     = data\$bkpBegin,
                    bkpEnd       = data\$bkpEnd,
                    bkpLength    = data\$bkpLength,
                    status       = data\$status)

dev.off()

write.table(output,
            file      = "${seg}",
            sep       = "\t",
            quote     = F,
            col.names = F,
            row.names = F,
            append    = F)

EndScript

  close(OUT);

  # Execute R with the script
  my $logR = "${scriptFile}.log.Rout";
  system("$R_CMD ${scriptFile} ${logR}");
  if (! -e $seg || ! -e $png) {
    print STDERR "Warning: Breakpoint number ${number} - Segmentation failed\n";
  } else {
    `mv ${seg} ${txtFile}`;
    if ($lastzSA ne "" || $lastzSB ne "") {
      `mv ${png} ${pngFile}`;
    } else {
      print STDERR "Warning: Breakpoint number ${number} - There are no hits (SR vs SA and SR vs SB)\n";
    }
  }
  deleteFile($scriptFile);
  deleteFile($logR);
  deleteFile($seg);
  deleteFile($png);
}
