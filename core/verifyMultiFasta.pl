#!/usr/bin/perl -w

###############################################################################
# Copyright (C) 2023 Inria                         #
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
require "${Bin}/util.pm";

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################

# Receive and process the parameters
my %PARAMETERS;
validateParameters();

# Read the FASTA directory and create the output file
readMultiFasta();

exit(0);

###############################################################################
## AUXILIARY SUBROUTINES ######################################################
###############################################################################

###############################################################################
# This function prints the script usage
sub printUsage {
  
  my ($error) = @_;

  if (defined $error) {
    print STDERR "\nERROR: $error\n\n";
  }

  print << "End_print_usage";

Usage:
  print verifyFASTAdirectory.pl <fasta> <outputfile>

Mandatory parameters:

 <fasta> Path to the multi-fasta file that contain all sequences of a genome.

 <outputfile> Name of the file that will receive the script result.

Script description:

  The script reads the multi-fasta file <fasta>
and produces a table which contains, for each FASTA entry, the following columns:
    + chr  - Name of the chromosome
    + inf  - Leftmost extremity of the chromosome (position 1)
    + sup  - Rightmost extremity of the chromosome (length of the sequence)
    + file - Name of the genome FASTA file


End_print_usage

  if (defined $error) {
    print "\nERROR: $error\n\n";
  }

  exit(1);
}

###############################################################################
# This function processes the parameters
sub validateParameters {

  my @parameters = @ARGV;
  my $nParameters = @parameters;
  
  if ($nParameters != 2) {
    printUsage("Invalid number of parameters.");
  }

  my $fasta  = fileExists(trim($parameters[0]), \&printUsage);
  my $outputFile = trim($parameters[1]);

  if ($outputFile eq "") {
    printUsage("Missing parameter <outputfile>");
  }
  $outputFile = getAbsolutePath($outputFile);
  
  $PARAMETERS{"FASTA"}  = $fasta;
  $PARAMETERS{"OUTPUTFILE"} = $outputFile;
}



###############################################################################
# This function reads the fasta file and returns the size of the sequence that
# is inside of it (if more than one sequence is found, the function returns -1)
sub readMultiFasta {

  my $fasta  = $PARAMETERS{"FASTA"};
  my $outputFile = $PARAMETERS{"OUTPUTFILE"};
  
  open(OUT, ">${outputFile}") or die printUsage("Could not write ${outputFile}");
  my $current_name = "";
  my $length = 0;
  my $nSequences = 0;
  open(IN, "${fasta}") or die printUsage("Could not read the ${fasta}");
  while (my $line = <IN>) {
    $line =~ s/\s+//g;
    if ($line !~ /^>/) {
      $length += length($line);
    } elsif ($line =~ /^>/) {
      if ($nSequences >0) {
      	# Treat previous chromosome :
      	print OUT "${current_name}\t1\t${length}\t${fasta}\n";
      }
      # update length and chromosome name with next sequence
      $length = 0;
      ($current_name) = $line =~ /^>(.+)$/;
      $current_name =~ tr/a-z/A-Z/;
      $nSequences++;
    }
  }
  close(IN);
  # writes last sequence
  print OUT "${current_name}\t1\t${length}\t${fasta}\n";
  
  close(OUT);

}
