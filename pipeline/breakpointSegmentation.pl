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
require "${Bin}/../core/util.pm";

###############################################################################
# Constants: important columns (chromosome and sequence boundaries)
my $BKP       = 0;
my $SR_CHR    = 4;
my $SR_INF    = 7;
my $BKP_BEGIN = 19;
my $BKP_END   = 20;
my $STATUS    = 21;

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################

# Verify the environment
my %ENVIRONMENT;
verifyEnvironment();

# Receive and process the parameters
my %PARAMETERS;
validateParameters();

# Read the fasta directory for getting the SR FASTA files
my %SR_FASTAS;
getFASTAFiles();

# Read the alignment directory for getting the LASTZ result files
my %LASTZ_FILES;
getLASTZFiles();

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
  perl breakpointSegmentation.pl <table> <fastadir> <aligndir> <outputdir> [tmpdir]

-------------------------------------------------------------------------------
Mandatory parameters:

  <table> Path for the file that contains the list of breakpoints.

  <fastadir> Path for the directory where the FASTA file of the sequences
             SR, SA and, SB can be found.

  <aligndir> Path for the directory where the script will write the results
             of the alignment of the sequences SR against SA and SR against
             SB.

  <outputdir> Path for the directory where the script will write the result 
              of the segmentation of each breakpoint.

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

-------------------------------------------------------------------------------
SCRIPT OUTPUT:

  The script will generate the segmentation information and write it in the
  directory <outputdir>. For each breakpoint, the script will create two 
  files:
  - Image file    : Breakpoint_[N].png
  - Data file     : Breakpoint_[N].txt
  where: [N] is the number of the breakpoint.

  The image file contains the graphical representation of the segmentation.
  The data file contains the information about the segmentation. It contains,
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
  
  if ($nParameters < 4 && $nParameters > 5) {
    printUsage("Invalid number of parameters.");
  }

  my $table     = fileExists(trim($parameters[0]), \&printUsage);
  my $fastadir  = directoryExists(trim($parameters[1]), \&printUsage);
  my $aligndir  = directoryExists(trim($parameters[2]), \&printUsage);
  my $outputdir = directoryExists(trim($parameters[3]), \&printUsage);
  my $tmpDir    = "/tmp";
  if ($nParameters == 5) {
    $tmpDir = directoryExists(trim($parameters[4]), \&printUsage);
  }

  $PARAMETERS{"TABLE"}     = $table;
  $PARAMETERS{"FASTADIR"}  = $fastadir;
  $PARAMETERS{"ALIGNDIR"}  = $aligndir;
  $PARAMETERS{"OUTPUTDIR"} = $outputdir;
  $PARAMETERS{"TMPDIR"}    = $tmpDir;

}

###############################################################################
# This function verifies if the core script is ok
sub verifyEnvironment {
  my $corescript = "${Bin}/../core/segmentation.pl";
  if (! -e $corescript) {
    printUsage("Missing script ${corescript}");
  }
  $ENVIRONMENT{"CORESCRIPT"} = $corescript;
}

###############################################################################
# This function get the list of FASTA files of the sequences SR
sub getFASTAFiles {

  my $fastaDir = $PARAMETERS{"FASTADIR"};
  my @files = `ls -1 ${fastaDir}/`;
  chomp(@files);
  
  foreach my $file (@files) {
    if ($file =~ /^Breakpoint_(\d+)_SR\.fasta$/) {
      my $number = $1;
      my $key    = sprintf("%04d", $number);
      $SR_FASTAS{$key} = "${fastaDir}/${file}";
    }
  }
  
}

###############################################################################
# This function get the list of FASTA files of the sequences SR
sub getLASTZFiles {

  my $alignDir = $PARAMETERS{"ALIGNDIR"};
  my @files = `ls -1 ${alignDir}/`;
  chomp(@files);

  my %discarded;
  
  foreach my $file (@files) {
    if ($file =~ /^Breakpoint_(\d+)_S[AB]\.lastz$/) {
      
      my $number = $1;
      my $key    = sprintf("%04d", $number);
      
      if (! defined $LASTZ_FILES{$key} && ! defined $discarded{$key}) {
	
	my $sa = "${alignDir}/Breakpoint_${number}_SA.lastz";
	my $sb = "${alignDir}/Breakpoint_${number}_SB.lastz";
	
	my $ok = 1;
	if (! -e $sa) {
	  print "Warning: Missing file ${sa} [discarding ${file}]\n";
	  $ok = 0;
	}
	if ($ok == 1 && ! -e $sb) {
	  print "Warning: Missing file ${sb} [discarding ${file}]\n";
	  $ok = 0;
	}

	if ($ok == 1) {
	  $LASTZ_FILES{$key} = "$sa;$sb";
	} else {
	  $discarded{$key} = 1;
	}

      }
      
    }

  }
  
}


###############################################################################
# This function calls the segmentation script for each breakpoint region
sub segmentation {

  my $table     = $PARAMETERS{"TABLE"};
  my $outputdir = $PARAMETERS{"OUTPUTDIR"};

  # Get the list of breakpoints
  my %statusValues;
  my %data;
  open(IN, "${table}") or die printUsage("Could not read the file ${table}");
  while (my $line = <IN>) {
    chomp($line);

    my @columns = split("\t", $line);
    my $id       = $columns[$BKP];
    my $key      = sprintf("%04d", $id);
    my $chr      = $columns[$SR_CHR];
    my $sRinf    = $columns[$SR_INF];
    my $bkpBegin = $columns[$BKP_BEGIN];
    my $bkpEnd   = $columns[$BKP_END];
    my $status   = $columns[$STATUS];
    
    my $length = $bkpEnd - $bkpBegin;
    my $begin  = $sRinf + $bkpBegin;
    my $end    = $sRinf + $bkpEnd;

    $statusValues{$key} = $status;
    $data{$key} = join("\t", ($id, $chr, $begin, $end, $length, 
			      $begin, $end, $length));    
  }
  close(IN);

  # Perform the segmentation
  foreach my $number (sort keys %LASTZ_FILES) {
    if (defined $SR_FASTAS{$number}) {
      my $fasta = $SR_FASTAS{$number};
      my ($sa, $sb) = split(";", $LASTZ_FILES{$number});
      print "Segmentation Breakpoint ${number}\n";
      runScript($number, $sa, $sb, $fasta);
    }
  }

  # Verify the segmentation result
  foreach my $number (sort keys %statusValues) {
    my $status = $statusValues{$number};
    my $breakpointData = $data{$number};
    my $txtFile = "${outputdir}/Breakpoint_${number}.txt";
    if ($status == 1 && ! -e $txtFile) {
      $status = -7;
    }
    if ($status != 1) {
      open(OUT, ">${txtFile}");
      print OUT "$breakpointData\t$status\n";
      close(OUT);
    }
  } 
}

###############################################################################
# This function runs the core script with the received parameters
sub runScript {

  my ($number, $sa, $sb, $fasta) = @_;

  my $script    = $ENVIRONMENT{"CORESCRIPT"};
  my $table     = $PARAMETERS{"TABLE"};
  my $outputdir = $PARAMETERS{"OUTPUTDIR"};
  my $tmpDir    = $PARAMETERS{"TMPDIR"};

  my $txtFile = "${outputdir}/Breakpoint_${number}.txt";
  my $pngFile = "${outputdir}/Breakpoint_${number}.png";

  my $cmd = "perl ${script} ${number} ${table} ${sa} ${sb} ${fasta} ${txtFile} ${pngFile} ${tmpDir}";
  system($cmd);

}
