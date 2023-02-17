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

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################

###############################################################################
# Process the received parameters
my ($inputfile, $outputfile) = processParameters();

###############################################################################
# Read the sequence and identify the masked intervals
my $sequence = "";
my $nSequences = 0;

my $accumulated = 0;
my $intervalBegin = 0;

open(IN, "$inputfile") or die printUsage("Could not read the file $inputfile.");
open(OUT, ">$outputfile") or die printUsage("Could not write the file $outputfile.");
while (my $line = <IN>) {
  $line =~ s/\s+//g;
  if ($line !~ /^>/ && $line ne "") {
      
    if ($line =~ /^[^nN]/) {
      $line =~ s/^([^nN]+)(.*)$/$2/;
      my $size = length($1);
      if ($intervalBegin != 0) {
	print OUT "$intervalBegin\t$accumulated\n";
	$intervalBegin = 0;
      }
      $accumulated += $size;
    }

    while ($line =~ /^([N]+)(.*)$/) {
      $line = $2;
      my $size = length($1);
      if ($intervalBegin == 0) {
	$intervalBegin = $accumulated + 1;
      }
      $accumulated += $size;
      if ($line =~ /^[^nN]/) {
	$line =~ s/^([^nN]+)(.*)$/$2/;
	$size = length($1);
	if ($intervalBegin != 0) {
	  print OUT "$intervalBegin\t$accumulated\n";
	  $intervalBegin = 0;
	}
	$accumulated += $size;
      }
    }

  } elsif ($line =~ /^>/) {
    $nSequences++;
    if ($nSequences > 1) {
      printUsage("$inputfile has more than one sequence.");
    }
  }

}

if ($intervalBegin != 0) {
  print OUT "$intervalBegin\t$accumulated\n";
}

close(IN);

exit(0);

###############################################################################
## AUXILIARY SUBROUTINES ######################################################
###############################################################################
# This function processes the parameters
sub processParameters {
  
  my ($inputFile, $outputFile) = @ARGV;

  if (defined $inputFile) {
    $inputFile = getAbsolutePath($inputFile);
    if (! -e $inputFile) {
      printUsage("Could not find the file $inputFile.");
    }
  } else {
    printUsage("Missing parameter <fastafile>");
  }

  if (! defined $outputFile) {
    printUsage("Missing parameter <outputfile>");
  }
  $outputFile = getAbsolutePath($outputFile);
  
  return ($inputFile, $outputFile);
}

###############################################################################
# This function prints the script usage
sub printUsage {

  my ($error) = @_;

  if (defined $error) {
    print STDERR "\nERROR: $error\n\n";
  }
  print << "End_print_usage";

Usage:
  print getMaskedIntervals.pl <fastafile> <outputfile>

Mandatory parameters:

 <fastafile> Fasta file that has the sequence that will be analyzed.

 <outputfile> File that will receive the script result.

Script description:

  The script reads the FASTA file and produces a table that have two columns:
begin and end coordinates of the masked intervals (intervals with N). The FASTA
file must have at most one sequence.

End_print_usage

  if (defined $error) {
    print "\nERROR: $error\n\n";
  }

  exit(1);
}

###############################################################################
# This function returns the absolute path of a file
sub getAbsolutePath {
  my ($relativePath) = @_;
  return File::Spec->rel2abs($relativePath);
}
