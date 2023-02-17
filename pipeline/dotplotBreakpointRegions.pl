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

# Generate the dotplots
dotPlotBreakpointRegions();

exit(0);

###############################################################################
## AUXILIARY SUB-ROUTINES #####################################################
###############################################################################
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
  perl dotplotBreakpointRegions.pl <table> <fastadir> <aligndir> <outputdir> [tmpdir]

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

  The script will generate a dotplot for each breakpoint and write them in the
  directory <outputdir>. The dotplot file has the following name pattern:
  - Image file    : Breakpoint_[N].png
  where: [N] is the number of the breakpoint.

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
  my $corescript = "${Bin}/../core/dotplotBreakpoint.pl";
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
# This function calls the script that generate the dotplot for each breakpoint
sub dotPlotBreakpointRegions {

  foreach my $number (sort keys %LASTZ_FILES) {
    if (defined $SR_FASTAS{$number}) {
      my $fasta = $SR_FASTAS{$number};
      my ($sa, $sb) = split(";", $LASTZ_FILES{$number});
      print "Dotplot Breakpoint ${number}\n";
      runScript($number, $sa, $sb, $fasta);
    }
  }

}

###############################################################################
# This function run the core script with the received parameters
sub runScript {

  my ($number, $sa, $sb, $fasta) = @_;

  my $script    = $ENVIRONMENT{"CORESCRIPT"};
  my $table     = $PARAMETERS{"TABLE"};
  my $outputdir = $PARAMETERS{"OUTPUTDIR"};
  my $tmpDir    = $PARAMETERS{"TMPDIR"};

  my $pngFile = "${outputdir}/Breakpoint_${number}.png";

  my $cmd = "perl ${script} ${number} ${table} ${sa} ${sb} ${fasta} ${pngFile} ${tmpDir}";
  system($cmd);

}
