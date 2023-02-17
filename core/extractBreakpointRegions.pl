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
# Constants: important columns (chromosome and sequence boundaries)
my $BKP    = 0;
my $SR_CHR = 4;
my $SR_INF = 7;
my $SR_SUP = 8;
my $SA_CHR = 11;
my $SB_CHR = 12;
my $SA_INF = 15;
my $SA_SUP = 16;
my $SB_INF = 17;
my $SB_SUP = 18;
my $STATUS = 21;

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################

# Receive and process the parameters
my %PARAMETERS;
validateParameters();

# Read the FASTA file
my $SEQUENCE = "";
readFASTA();

# Read the table and create the breakpoint sequences
createBreakpointFASTAs();

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
  perl extractBreakpointRegions.pl <table> <chromosome> <fasta> <outputdir> <genome>

-------------------------------------------------------------------------------
Mandatory parameters:

  <table> Path for the file that contains the list of breakpoints.

  <chromosome> Name of the chromosome that is related to the given FASTA file.

  <fasta> Path for the file that has the DNA sequence of the chromosome.

  <outputdir> Path for the directory where the script will write the FASTA
              files of the breakpoint regions.

  <genome> Indicate the genome which owns the specified chromosome. This 
           parameter must be set with 1 to indicate the genome GR and with 2
           to indicate the genome GO.

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

  + FASTA file

  The FASTA file must have a single sequence that represents the whole
  DNA sequence of the chromosome.

-------------------------------------------------------------------------------
Script output:

The script will write inside of the directory <outputdir> all breakpoint 
sequences that are present in the chromosome <chromosome> and are listed on
the file <table>. Note that if the user wants to create the sequences SR, it 
must set the parameter <fasta> with the FASTA file of the desired chromosome 
of the genome GR and set the parameter <genome> with the value 1. If the user
wants to create the sequences SA and SB, it must set the parameter <fasta> 
with the FASTA file of the desired chromosome of the genome GO and set the 
parameter <genome> with the value 2.

The script will create FASTA files that have the following name pattern:
Breakpoint_[N]_S[X].fasta, where [N] is the number of the breakpoint and [X]
is one of the letters R, A or, B.

End_print_usage

  if (defined $error) {
    print "\nERROR: $error\n\n";
  }
  
  exit(1);
}

###############################################################################
# This function validate the received parameters
sub validateParameters {
  
  my @parameters = @ARGV;
  my $nParameters = @parameters;
  
  if ($nParameters != 5) {
    printUsage("Invalid number of parameters.");
  }
  
  my $table     = fileExists(trim($parameters[0]), \&printUsage);
  my $chr       = trim($parameters[1]);
  my $fasta     = fileExists(trim($parameters[2]), \&printUsage);
  my $outputdir = directoryExists(trim($parameters[3]), \&printUsage);
  my $genome    = trim($parameters[4]);
  
  $chr =~ tr/a-z/A-Z/;
  
  if ($genome !~ /^[12]$/) {
    printUsage("Invalid value for the parameter genome. It must be set with 1 or 2 (GR or GO).");
  }
  
  $PARAMETERS{"TABLE"}     = $table;
  $PARAMETERS{"CHR"}       = $chr;
  $PARAMETERS{"FASTA"}     = $fasta;
  $PARAMETERS{"OUTPUTDIR"} = $outputdir;
  $PARAMETERS{"GENOME"}    = $genome;
}

###############################################################################
# This function reads the FASTA file
sub readFASTA {
  my $fastaFile = $PARAMETERS{"FASTA"};  
  my $nSequences = 0;
  open(IN, "${fastaFile}") or die printUsage("Could not read the file ${fastaFile}");
  while (my $line = <IN>) {
    $line =~ s/\s+//g;
    if ($line !~ /^>/ && $nSequences == 1) {
      $SEQUENCE .= $line;
    } elsif ($line =~ /^>/) {
      $nSequences++;
      if ($nSequences > 1) {
	printUsage("${fastaFile} has more than one sequence.");
      }
    }
  }
  close(IN);
}

###############################################################################
# This function reads the breakpoint table and creates the FASTA files that
# of fragments that are in the specified chromosome.
sub createBreakpointFASTAs {

  my $table     = $PARAMETERS{"TABLE"};
  my $chr       = $PARAMETERS{"CHR"};
  my $outputDir = $PARAMETERS{"OUTPUTDIR"};
  my $genome    = $PARAMETERS{"GENOME"};

  open(IN, "${table}") or die printUsage("Could not read the file ${table}");
  while (my $line = <IN>) {

    chomp($line);
    my @columns = split("\t", $line);

    my $breakpointNumber = sprintf("%04d", $columns[$BKP]);
    my $status = $columns[$STATUS];
    if ($status == 1) {
    
      if ($genome == 1) {
	# Genome GR - Fragment SR
	my $chrGR = $columns[$SR_CHR];
	$chrGR =~ tr/a-z/A-Z/;
      
	if ($chr eq $chrGR) {
	  my $name = "Breakpoint_${breakpointNumber}_SR";
	  my $fileName = "${outputDir}/${name}.fasta";
	  writeFASTA($name, $fileName, $columns[$SR_INF], $columns[$SR_SUP]);
	}

      } else {
	# Genome GO - Fragments SA and SB
      
	# SA
	my $chrA = $columns[$SA_CHR];
	$chrA =~ tr/a-z/A-Z/;
	if ($chr eq $chrA) {
	  my $name = "Breakpoint_${breakpointNumber}_SA";
	  my $fileName = "${outputDir}/${name}.fasta";
	  writeFASTA($name, $fileName, $columns[$SA_INF], $columns[$SA_SUP]);
	}

	# SB
	my $chrB = $columns[$SB_CHR];
	$chrB =~ tr/a-z/A-Z/;
	if ($chr eq $chrB) {
	  my $name = "Breakpoint_${breakpointNumber}_SB";
	  my $fileName = "${outputDir}/${name}.fasta";
	  writeFASTA($name, $fileName, $columns[$SB_INF], $columns[$SB_SUP]);
	}
      }
    }
  }
  close(IN);

}

###############################################################################
# This function writes the FASTA file
sub writeFASTA {
  my ($name, $file, $begin, $end) = @_;
  my $length = $end - $begin + 1;
  my $fragment = substr($SEQUENCE, $begin - 1, $length);
  open(OUT, ">${file}") or die printUsage("Could not create the FASTA file ${file}");
  print OUT ">$name\n";
  my $index = 0;
  while ($index < $length) {
    my $size = 10000;
    if ($index + $size > $length) {
      $size = $length - $index;
    }
    my $substring = substr($fragment, $index, $size);
    print OUT "$substring\n";
    $index += $size;
  }
  close(OUT);
}
