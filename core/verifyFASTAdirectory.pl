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

# Read the FASTA directory and create the output file
readFASTADirectory();

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
  print verifyFASTAdirectory.pl <directory> <outputfile>

Mandatory parameters:

 <directory> Name of the directory where the script can find the FASTA files.

 <outputfile> Name of the file that will receive the script result.

Script description:

  The script read all FASTA files that are located in the directory <directory>
and produces a table which contains, for each FASTA file, the following columns:
    + chr  - Name of the chromosome
    + inf  - Leftmost extremity of the chromosome (position 1)
    + sup  - Rightmost extremity of the chromosome (lenght of the sequence)
    + file - Name of the chromosome FASTA file

  Warning: The FASTA file must have just one sequence. Otherwise, it will be
           ignored by the script.

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

  my $directory  = directoryExists(trim($parameters[0]), \&printUsage);
  my $outputFile = trim($parameters[1]);

  if ($outputFile eq "") {
    printUsage("Missing parameter <outputfile>");
  }
  $outputFile = getAbsolutePath($outputFile);
  
  $PARAMETERS{"DIRECTORY"}  = $directory;
  $PARAMETERS{"OUTPUTFILE"} = $outputFile;
}


###############################################################################
# This function reads the FASTA directory
sub readFASTADirectory {

  my $directory  = $PARAMETERS{"DIRECTORY"};
  my $outputFile = $PARAMETERS{"OUTPUTFILE"};
  
  # Get the list of files of the directory
  my @files = `ls -1 ${directory}/`;
  chomp(@files);

  open(OUT, ">${outputFile}") or die printUsage("Could not write ${outputFile}");
  foreach my $file (sort @files) {
    if ($file =~ /^(.+)\.fasta$/i) {
      my $chr = $1;
      $chr =~ tr/a-z/A-Z/;
      my $length = readFASTA("${directory}/${file}");
      if ($length > 0) {
	print OUT "${chr}\t1\t${length}\t${directory}/${file}\n";
      } elsif ($length == 0) {
	print "Warning: ${directory}/${file} has no sequence.\n";
      } else {
	print "Warning: ${directory}/${file} has more than one sequence.\n";
      }
    }
  }
}
close(OUT);

###############################################################################
# This function reads the fasta file and returns the size of the sequence that
# is inside of it (if more than one sequence is found, the function returns -1)
sub readFASTA {
  my ($file) = @_;
  my $length = 0;
  my $nSequences = 0;
  open(IN, "${file}") or die printUsage("Could not read the ${file}");
  while (my $line = <IN>) {
    $line =~ s/\s+//g;
    if ($line !~ /^>/ && $nSequences == 1) {
      $length += length($line);
    } elsif ($line =~ /^>/) {
      $nSequences++;
      if ($nSequences > 1) {
	return -1;
      }
    }
  }
  close(IN);
  return $length;
}
