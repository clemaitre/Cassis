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
# Get the parameters from the command line
my ($orthologousTable, $secondSpeciesGeneTable,
    $chromosomesToIgnore1, $chromosomesToIgnore2,
    $includeApparent) = @ARGV;

# Verify the parameters
verifyParameters();

###############################################################################
# Arrays which will keep the list of regular expressions to discard chromosomes
my @toIgnore1;
my @toIgnore2;

# Hash that will keep the map between gene name and strand
my %geneStrand;

###############################################################################
## Main Program ###############################################################
###############################################################################

# Read the table of genes of the second species
open(IN, "$secondSpeciesGeneTable")
  or die "Could not open the file $secondSpeciesGeneTable\n";

while (my $line = <IN>) {
  chomp($line);
  if ($line !~ /^Ensembl/) {
    my ($id, $chr, $start, $end, $strand) = split("\t", $line);
    $geneStrand{$id} = $strand;
  }
}

close(IN);

# Read the ortology, add the information about strand of the second species
# and output (to STDOUT)
open(IN, "$orthologousTable")
  or die "Could not open the file $orthologousTable\n";

my %repeated1;
my %repeated2;

# Print the header
while (my $line = <IN>) {

  chomp($line);
  # Ignore the header of the file
  if ($line !~ /^Ensembl/) {

    my ($id1, $chr1, $start1, $end1, $strand1, $id2, $chr2,
	$start2, $end2, $orthology) = split("\t", $line);

    # Consider just the one2one genes
    my $ok = 0;
    if ($includeApparent == 1) {
      if ($orthology =~ /one2one/i) {
	$ok = 1;
      }
    } else {
      if ($orthology =~ /^ortholog_one2one$/i) {
	$ok = 1;
      }
    }

    if ($ok == 1) {
      if (defined $geneStrand{$id2}) {
	if (discard($chr1, @toIgnore1) == 0 && discard($chr2, @toIgnore2) == 0) {
	  if (! defined $repeated1{$id1} && ! defined $repeated1{$id2}) {
	    my @out = ($id1, $chr1, $start1, $end1, $strand1, $id2, $chr2,
		       $start2, $end2, $geneStrand{$id2});
	    print join("\t", @out)."\n";
	    $repeated1{$id1} = 1;
	    $repeated2{$id2} = 1;
	  } else {
	    if (defined $repeated1{$id1}) {
	      print STDERR "WARNING: The gene $id1 from the first species appears more than one time.\n$line\n";
	    }
	    if (defined $repeated2{$id2}) {
	      print STDERR "WARNING: The gene $id2 from the second species appears more than one time.\n$line\n";
	    }
	  }
	}
      } else {
	print STDERR "$line\n";
	print STDERR "Error: $id2 does not appears in the file $secondSpeciesGeneTable\n";
      }
    }
  }
}

close(IN);


###############################################################################
## Auxiliary sub routines #####################################################
###############################################################################

###############################################################################
## Verify the parameters
sub verifyParameters {

  if (! defined $orthologousTable) {
    printUsage("Missing parameter <ortho_table>");
  }

  if (! defined $secondSpeciesGeneTable) {
    printUsage("Missing parameter <gene_table>");
  }

  $orthologousTable = fileExists($orthologousTable);

  $secondSpeciesGeneTable = fileExists($secondSpeciesGeneTable);

  if (! defined $chromosomesToIgnore1) {
    printUsage("Missing parameter <ignore1>");
  }

  if (! defined $chromosomesToIgnore2) {
    printUsage("Missing parameter <ignore2>");
  }

  if ($chromosomesToIgnore1 ne "") {
    @toIgnore1 = split(",", $chromosomesToIgnore1);
  }

  if ($chromosomesToIgnore2 ne "") {
    @toIgnore2 = split(",", $chromosomesToIgnore2);
  }

  if (! defined $includeApparent) {
    $includeApparent = 0;
  } else {
    if ($includeApparent ne "1" && $includeApparent ne "0") {
      print STDERR "\nERROR: Invalid value for the parameter <apparent>\n";
      printUsage();
    }
  }

}

###############################################################################
## Return 1 with the value $chr match at least one of the given
## regular expressions
sub discard {
  my ($chr, @regularExpressions) = @_;
  foreach my $re (@regularExpressions) {
    if ($chr =~ /$re/i) {
      return 1;
    }
  }
  return 0;
}

###############################################################################
## Print information about how to use the script
sub printUsage {

  my ($error) = @_;

  if (defined $error) {
    print STDERR "\nERROR: $error\n\n";
  }

print <<"End_Print_Usage";

-------------------------------------------------------------------------------
Usage:
perl createOne2OneTable.pl <ortho_table> <gene_table> <ignore1> <ignore2> <apparent>

PARAMETERS:

<ortho_table> Table that contains the orthology relationship between 2 species.

<gene_table>  Table that contanis the list of genes of the second species.

<ignore1>     List of regular expression (separated by comma ',') which will be
              used to discard undesired chromosomes from the first species.
              If no chromosome should be discarded, use empty string ("")

<ignore1>     List of regular expression (separated by comma ',') which will be
              used to discard undesired chromosomes from the second species.

<apparent>    [1/0] - If apparent = 1, the script will also preserve the pairs
              of ortholog genes that have the relationship of the type
              apparent_ortholog_one2one, if apparent = 0, just the pairs which
              are classified as ortholog_one2one will be preserved. [Default
              value = 0]

EXAMPLE:
perl createOne2OneTable.pl Hsapiens_Mmusculus.txt MmusculusGenes.txt "Y,MT" "" 0

In this example, the script will unify the given tables, selecting just the
genes that have an one2one ortholog relationship and discarding all genes have
the chromosome Y or MT in the first species (Homo sapiens in this case).

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
TABLE FORMATS:

The tables used by this program should be produced by Biomart (Ensembl).

-------------------------------------------------------------------------------
The table <ortho_table> contains the orthology relationship between 2 species.

To produce this table the user must access, on Biomart, the database of the
first species and them select Attributes -> Homologs and check the following
fields:

For the first species (Section GENE)
- Ensembl Gene ID
- Chromosome Name
- Gene Start (bp)
- Gene End (bp)
- Strand

For the second species (Section "Species" ORTHOLOGS - Ex: RAT ORTHOLOGS)
- "Species" Ensembl Gene ID
- "Species" Chromosome Name
- "Species" Gene Start (bp)
- "Species" Gene End (bp)
- Orthology Type

The columns of the table must respect this order. The user must perform the
search and save the file as a TSV -- TAB SEPARATED VALUES

-------------------------------------------------------------------------------
The table <ortho_table> contains the list of all genes of the second species.

To produce this table, the user must access, on Biomart, the database of the
second species and them select Attributes -> Features and check the following
fields:

Section GENE:
- Ensembl Gene ID
- Chromosome Name
- Gene Start (bp)
- Gene End (bp)
- Strand

The columns of the table must respect this order. The user must perform the
search and save the file as a TSV -- TAB SEPARATED VALUES

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
OUTPUT:

The script will output to the STDOUT a table that has all ortholog genes that
have relationship of the type one2one and that are in chromosomes that do not
appear in the list of chromosomes that should be ignored. The script also
ignores the cases where a gene appears on more the one relationship one2one
(it keeps just the first reference). The table has the following fields:

- g1      = Name of the gene in the first species
- c1      = Name of the chromosome in the first species
- inf1    = Start position of the gene in the first species
- sup1    = End position of the gene in the first species
- strand1 = Strand where the gene is located in the first species
- g2      = Name of the gene in the second species
- c2      = Name of the chromosome in the second species
- inf2    = Start position of the gene in the second species
- sup2    = End position of the gene in the second species
- strand2 = Strand where the gene is located in the second species

End_Print_Usage

  if (defined $error) {
    print "\nERROR: $error\n\n";
  }

  exit(1);
}

###############################################################################
# This function verifies if the file exists and returns its full path
sub fileExists {
  my ($relativePath) = @_;
  my $file = getAbsolutePath($relativePath);
  if (! -e $file) {
    printUsage("Could not find the file $file");
  }
  return $file;
}

###############################################################################
# This function returns the absolute path of a file
sub getAbsolutePath {
  my ($relativePath) = @_;
  return File::Spec->rel2abs($relativePath);
}
