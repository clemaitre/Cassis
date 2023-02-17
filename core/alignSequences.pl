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
## MAIN PROGRAM ###############################################################
###############################################################################

# Receive and process the parameters
my %PARAMETERS;
validateParameters();

# Align the sequences
runLASTZ();

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
  perl alignSequences.pl <fasta1> <fasta2> <outputfile> <lastzparameters> <matrix> <lastz> [tmpdir]

-------------------------------------------------------------------------------
Mandatory parameters:

  <fasta1> Path for the FASTA file that has the sequence SR.

  <fasta2> Path for the FASTA file that has the sequence SA (or SB).

  <outputfile> name of the file that will receive the results of the
               alignment.

  <lastzparameters> string that have the list o LASTZ parameters (Except
                    parameter Q)
                    Example: "K=3000 E=30 O=400 L=3000 H=2000 B=2"

  <matrix> matrix that will be used by the LASTZ to compute the alignment.
           If the LASTZ must use the default scoring values, use the value
           "default" for this parameter.

  <lastz> Path for the program LASTZ.

-------------------------------------------------------------------------------
Optional parameter:

  [tmpdir] Path for the temporary directory. /tmp/ is the default temporary 
           directory.


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
  
  my $fasta1     = fileExists(trim($parameters[0]), \&printUsage);
  my $fasta2     = fileExists(trim($parameters[1]), \&printUsage);
  my $output     = trim($parameters[2]);
  my $parameters = trim($parameters[3]);
  my $matrix     = trim($parameters[4]);
  my $lastz     = fileExists(trim($parameters[5]), \&printUsage);
  my $tmpDir   = "/tmp";
  if ($nParameters == 7) {
    $tmpDir = directoryExists(trim($parameters[6]), \&printUsage);
  }

  if ($matrix !~ /^default$/i) {
    $matrix = fileExists($matrix, \&printUsage);
  }
  
  $PARAMETERS{"FASTA1"}     = $fasta1;
  $PARAMETERS{"FASTA2"}     = $fasta2;
  $PARAMETERS{"OUTPUT"}     = $output;
  $PARAMETERS{"PARAMETERS"} = $parameters;
  $PARAMETERS{"MATRIX"}     = $matrix;
  $PARAMETERS{"LASTZ"}      = $lastz;
  $PARAMETERS{"TMPDIR"}     = $tmpDir;

}


###############################################################################
# This function runs the program lastz
sub runLASTZ {

  my $fasta1     = $PARAMETERS{"FASTA1"};
  my $fasta2     = $PARAMETERS{"FASTA2"};
  my $output     = $PARAMETERS{"OUTPUT"};
  my $parameters = $PARAMETERS{"PARAMETERS"};
  my $matrix     = $PARAMETERS{"MATRIX"};
  my $lastz      = $PARAMETERS{"LASTZ"};
  my $tmpDir     = $PARAMETERS{"TMPDIR"};

  my $tmpOutput = "${tmpDir}/".rand().time();

  my $command = "${lastz} ${fasta1} ${fasta2} ${parameters}";
  if ($matrix !~ /^default$/i) {
    $command .= " Q=${matrix} ";
  }
  $command .= " > ${tmpOutput}";
  system("$command");

  processBZOUT($tmpOutput, $output);
  deleteFile($tmpOutput);
}

###############################################################################
# This function process the lastz output
sub processBZOUT {
  
  my ($bzout, $lastzOutput) = @_;
  
  my $reversed = 0;
  my $length2nd = 0;
  my ($score, $b1, $b2, $e1, $e2);
  
  open(IN, $bzout) or die printUsage("Could not process the lastz output");
  open(OUT, ">${lastzOutput}") or die print("Could not write the output file");
  while (my $line = <IN>) {
    chomp($line);
    if ($line =~ /^s \{/) {
      <IN>;			#Jump the first sequence
      my $secondSequence = <IN>;
      chomp($secondSequence);
      if ($secondSequence =~ /^.+\-" \d+ (\d+) \d+ \d+$/) {
	$reversed = 1;
	$length2nd = $1;
      }
    } elsif ($line =~ /^\s+s\s+(\d+)/) {
      $score = $1;
    } elsif ($line =~ /^\s+b\s+(\d+)\s+(\d+)/) {
      $b1 = $1;
      if ($reversed == 0) {
	$b2 = $2;
      } else {
	$b2 = $length2nd - $2;
      }
    } elsif ($line =~ /^\s+e\s+(\d+)\s+(\d+)/) {
      $e1 = $1;
      if ($reversed == 0) {
	$e2 = $2;
      } else {
	$e2 = $length2nd - $2;
      }
      print OUT "$score\t$b1\t$b2\t$e1\t$e2\n";
    }
  }
  close(IN);
  close(OUT);
}
