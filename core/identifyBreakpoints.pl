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

# Run the script R
runScriptR();

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
  perl identifyBreakpoints.pl -i inputFile -t inputType -o outputFile [optional]

-------------------------------------------------------------------------------
Mandatory parameters:

  -i inputFile

     "inputFile" is the path for the file that has all pairs of orthologous
     genes  or synteny blocks between two genomes GR and GO. Note that only
     one2one orthologous pairs are allowed (see more information in the section
     "File formats").

  -t inputType

     "inputType" is the type of the input table: G for orthologous genes and
     B for synteny blocks.

  -o outputFile

     "outputFile" is the path for the file that will be written with the
     results of this script (see  more information in the section "File
     formats").

-------------------------------------------------------------------------------
Optional parameters:

  -R chromosomeBoundariesFileGR

     "chromosomeBoudariesFileGR" is the path for the file that has the
     begin and end position of all chromosomes of the genome GR. This
     data will be used to correct the extremities of the sequence R when
     the extension of the sequence falls outside of the chromosome
     sequence (see more information in the section "File formats").

  -O chromosomeBoundariesFileGO

     "chromosomeBoudariesFileGO" is the path for the file that has the
     begin and end position of all chromosomes of the genome GO. This
     data will be used to correct the extremities of the sequences A and B
     when they are related to a gene that is in the begin or in the end of
     a chromosome of the genome GO (see more information in the section
     "File formats").

  -C centromereBoundariesFile

     "centromereBoundariesFile" is the path for the file that has the begin
     and end position of the centromeres of all chromosomes of the genome
     GR. This data will be used to eliminate breakpoint regions that
     have a centromere inside of them (see more information in the section
     "File formats").

  -e [T/F]

     Extend - This parameter determine if the sequences SR, SA and, SB must
     include (T) or not (F) the genes A and B. This option is valid only for
     orthologous genes table (-i G)
     [Default value = T]

  -E [T/F]

     Extend previous/after - This parameter determine if the sequences SA
     and SB must include (T) or not (F) the genes that are located
     before/after the genes A/B. This option is valid only for orthologous
     genes table (-i G)
     [Default value = T]

  -b N

     Extend - This parameter determine if the sequences SR, SA and, SB must
     include a fragment of size N of the blocks that define the sequences.
     This option is valid only for synteny blocks table (-i B)
     [Default value = 50000]

  -B N

     Extend previous/after - This parameter determine if the sequences SA
     and SB must include a fragment of size N of the blocks that come
     before/after the blocks that define the sequences. This option is valid
     only for synteny blocks table (-i B)
     [Default value = 50000]

  -M N

     Minimum sequence size - N is a positive integer that indicates the
     minimum size that the sequences SR, SA and SB must have.
     [Default value for genes = 1]
     [Default value for synteny blocks = 50000]

  -v [T/F]

     Extend before verify length - This parameter determine if we will
     verify the minimum sequence size after (T) or before (F) applying
     the sequence extension.
     [Default value = T]

  -l N

     Maximum length of the sequence R - N is a integer that indicates the
     maximum number of bases of the sequence R.
     [Default value = 1000000000]

  -L N

     Maximum length of the sequences A/B - N is a integer that indicates
     the maximum number of bases of the sequences A/B.
     [Default value = 3000000]

  -n nonOverlappingIntervalsFile

     Name of the file that will receive the list of genes and merged intervals
     that was produced after resolving the problem of gene overlapping and was
     used to define the list of synteny blocks. This option is valid only for
     orthologous genes table (-i G)

  -m mergedIntervalsComposition

     Name of the file that will describe for each merged interval, the list of
     genes that is inside of it. This option is valid only for orthologous
     genes table (-i G)

  -s syntenyBlocksFile

     Name of the file that will receive the list of synteny blocks that were
     used to calculate the breakpoint regions. This option is valid only for
     orthologous genes table (-i G)

  -T temporaryDirectory

     Name of the temporary directory. /tmp/ is the default temporary directory.

-------------------------------------------------------------------------------
File formats:

  + Input File:

  Table of orthologous pairs between two genomes GR and GO. This file must
  have the following columns:

      + g1      - Name of the gene/block in the first genome (GR)
      + c1      - Chromosome where the gene/block g1 is located
      + inf1    - Start position of the gene/block g1 in the chromosome c1
      + sup1    - End position of the gene/block g1 in the chromosome c1
      + strand1 - Strand of the gene/block g1
      + g2      - Name of the gene/block in the second genome  (GO)
      + c2      - Chromosome where the gene/block g2 is located
      + inf2    - Start position of the gene/block g2 in the chromosome c2
      + sup2    - End position of the gene/block g2 in the chromosome c2
      + strand2 - Strand of the gene/block g2

  Warning: The file must NOT have header with the column names.

  + Chromosome Boundaries File:

  Table that have the information about the start and end position of each
  chromosome of the genome. It must have the following columns:

      + chr  - Name of the chromosome
      + inf  - Leftmost extremity of the chromosome
      + sup  - Rightmost extremity of the chromosome
      + file - Name of the chromosome FASTA file

  Warning: The file must NOT have header with the column names.

  + Centromere Boundaries File:

  Table that have the information about the start and end position of the
  centromere of each chromosome of the genome GR. This file must have
  the following columns:

      + chr - Name of the chromosome
      + inf - Leftmost extremity of the centromere of the chromosome
      + sup - Rightmost extremity of the centromere of the chromosome

  Warning: The file must NOT have header with the column names.

  + Output File:

  The script will output a table that contains information about the breakpoint
  regions that were indentified. This file has the following columns:

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

End_print_usage

  if (defined $error) {
    print "\nERROR: $error\n\n";
  }

  exit(1);

}

###############################################################################
# This function processes the parameters
sub validateParameters {

  my $nParameters = @ARGV;

  if ($nParameters % 2 != 0 || $nParameters == 0) {
    printUsage("Invalid number of parameters");
  }

  # Put the parameters on a hashtable
  for (my $i = 0; $i < $nParameters; $i+=2) {
    my $name = trim($ARGV[$i]);
    my $value = trim($ARGV[$i + 1]);
    if ($name !~ /^\-[ioeElLCmnsTtbBROMv]$/) {
      printUsage("Unknown parameter ($name)");
    }
    $PARAMETERS{$name} = $value;
  }

  if (! defined $PARAMETERS{"-i"}) {
    printUsage("Missing parameter -i (input file)");
  } else {
    $PARAMETERS{"-i"} = fileExists($PARAMETERS{"-i"}, \&printUsage);
  }

  if (! defined $PARAMETERS{"-t"}) {
    printUsage("Missing parameter -t (input type)");
  } else {
    my $type = $PARAMETERS{"-t"};
    if ($type !~ /^[gGbB]$/) {
      printUsage("Invalid value for the parameter <type>");
    } else {
      $type =~ tr/gb/GB/;
    }
    $PARAMETERS{"-t"} = $type;
  }

  if (! defined $PARAMETERS{"-o"}) {
    printUsage("Missing parameter -o (output file)");
  }
  $PARAMETERS{"-o"} = getAbsolutePath($PARAMETERS{"-o"});


  if (! defined $PARAMETERS{"-R"}) {
    $PARAMETERS{"-R"} = "";
  } else {
    $PARAMETERS{"-R"} = fileExists($PARAMETERS{"-R"}, \&printUsage);
  }

  if (! defined $PARAMETERS{"-O"}) {
    $PARAMETERS{"-O"} = "";
  } else {
    $PARAMETERS{"-O"} = fileExists($PARAMETERS{"-O"}, \&printUsage);
  }

  if (! defined $PARAMETERS{"-C"}) {
    $PARAMETERS{"-C"} = "";
  } else {
    $PARAMETERS{"-C"} = fileExists($PARAMETERS{"-C"}, \&printUsage);
  }

  if (! defined $PARAMETERS{"-e"}) {
    $PARAMETERS{"-e"} = "T";
  } else {
    my $value = $PARAMETERS{"-e"};
    $value =~ tr/a-z/A-Z/;
    if ($value ne "T" && $value ne "F") {
      printUsage("Invalid value for the parameter -e");
    }
    $PARAMETERS{"-e"} = $value;
  }

  if (! defined $PARAMETERS{"-E"}) {
    $PARAMETERS{"-E"} = "T";
  } else {
    my $value = $PARAMETERS{"-E"};
    $value =~ tr/a-z/A-Z/;
    if ($value ne "T" && $value ne "F") {
      printUsage("Invalid value for the parameter -E");
    }
    $PARAMETERS{"-E"} = $value;
  }

  if (! defined $PARAMETERS{"-b"}) {
    $PARAMETERS{"-b"} = 50000;
  } else {
    my $value = $PARAMETERS{"-b"};
    if ($value !~ /^\d+$/) {
      printUsage("Invalid value for the parameter -b");
    }
  }

  if (! defined $PARAMETERS{"-B"}) {
    $PARAMETERS{"-B"} = 50000;
  } else {
    my $value = $PARAMETERS{"-B"};
    if ($value !~ /^\d+$/) {
      printUsage("Invalid value for the parameter -B");
    }
  }

  if (! defined $PARAMETERS{"-M"}) {
    if ($PARAMETERS{"-t"} eq "B") {
      $PARAMETERS{"-M"} = 50000;
    } else {
      $PARAMETERS{"-M"} = 1;
    }
  } else {
    my $value = $PARAMETERS{"-M"};
    if ($value !~ /^\d+$/) {
      printUsage("Invalid value for the parameter -M");
    }
  }

  if (! defined $PARAMETERS{"-v"}) {
    $PARAMETERS{"-v"} = "T";
  } else {
    my $value = $PARAMETERS{"-v"};
    $value =~ tr/a-z/A-Z/;
    if ($value ne "T" && $value ne "F") {
      printUsage("Invalid value for the parameter -v");
    }
    $PARAMETERS{"-v"} = $value;
  }

  if (! defined $PARAMETERS{"-l"}) {
    $PARAMETERS{"-l"} = 1000000000;
  } else {
    my $value = $PARAMETERS{"-l"};
    $value =~ tr/A-Z/a-z/;
    if ($value !~ /^\d+$/ && $value !~ /^\d+e\d+$/) {
      printUsage("Invalid value for the parameter -l");
    }
  }

  if (! defined $PARAMETERS{"-L"}) {
    $PARAMETERS{"-L"} = 3000000;
  } else {
    my $value = $PARAMETERS{"-L"};
    $value =~ tr/A-Z/a-z/;
    if ($value !~ /^\d+$/ && $value !~ /^\d+e\d+$/) {
      printUsage("Invalid value for the parameter -L");
    }
  }

  if (! defined $PARAMETERS{"-m"}) {
    $PARAMETERS{"-m"} = "";
  } else {
    $PARAMETERS{"-m"} = getAbsolutePath($PARAMETERS{"-m"});
  }

  if (! defined $PARAMETERS{"-n"}) {
    $PARAMETERS{"-n"} = "";
  } else {
    $PARAMETERS{"-n"} = getAbsolutePath($PARAMETERS{"-n"});
  }

  if (! defined $PARAMETERS{"-s"}) {
    $PARAMETERS{"-s"} = "";
  } else {
    $PARAMETERS{"-s"} = getAbsolutePath($PARAMETERS{"-s"});
  }

  if (! defined $PARAMETERS{"-T"}) {
    $PARAMETERS{"-T"} = "/tmp";
  } else {
    $PARAMETERS{"-T"} = directoryExists($PARAMETERS{"-T"}, \&printUsage);
  }

}

###############################################################################
# This function verifies with the R scripts are ok
sub verifyEnvironment {

  my $overlappingScript = fileExists("${Bin}/R/overlappingGenes.R", \&printUsage);
  my $syntenyK2Script   = fileExists("${Bin}/R/synteny-k2-one2one.R", \&printUsage);
  my $breakpointsScript = fileExists("${Bin}/R/breakpoints.R", \&printUsage);
  my $utilScript        = fileExists("${Bin}/R/util.R", \&printUsage);
  my $intervalsScript   = fileExists("${Bin}/R/intervals.R", \&printUsage);
  $ENVIRONMENT{"OVERLAPPING"} = $overlappingScript;
  $ENVIRONMENT{"SYNTENY"}     = $syntenyK2Script;
  $ENVIRONMENT{"BREAKPOINTS"} = $breakpointsScript;
  $ENVIRONMENT{"UTIL"}        = $utilScript;
  $ENVIRONMENT{"INTERVAL"}    = $intervalsScript;

}

###############################################################################
# This function runs the script R with the received parameters
sub runScriptR {

  my $tableFile                   = $PARAMETERS{"-i"};
  my $type                        = $PARAMETERS{"-t"};
  my $outputFile                  = $PARAMETERS{"-o"};
  my $chromosomeFileGR            = $PARAMETERS{"-R"};
  my $chromosomeFileGO            = $PARAMETERS{"-O"};
  my $centromereFile              = $PARAMETERS{"-C"};
  my $mergedIntervalsComposition  = $PARAMETERS{"-m"};
  my $nonOverlappingIntervalsFile = $PARAMETERS{"-n"};
  my $syntenyBlocksFile           = $PARAMETERS{"-s"};
  my $limitR                      = $PARAMETERS{"-l"};
  my $limitAB                     = $PARAMETERS{"-L"};
  my $tempDir                     = $PARAMETERS{"-T"};
  my $minimumSequenceSize         = $PARAMETERS{"-M"};
  my $extendBeforeVerifyLength    = $PARAMETERS{"-v"};

  my $extended                    = $PARAMETERS{"-e"};
  my $extendedAB                  = $PARAMETERS{"-E"};
  if ($type eq "B") {
    $extended                     = $PARAMETERS{"-b"};
    $extendedAB                   = $PARAMETERS{"-B"};
  }

  my $overlappingScript = $ENVIRONMENT{"OVERLAPPING"};
  my $syntenyK2Script   = $ENVIRONMENT{"SYNTENY"};
  my $breakpointsScript = $ENVIRONMENT{"BREAKPOINTS"};
  my $utilScript        = $ENVIRONMENT{"UTIL"};
  my $intervalsScript   = $ENVIRONMENT{"INTERVAL"};

  my $scriptFile = $tempDir."/".rand().time();

  my $mergedIntervalsFile = "logFileName = NULL";
  if ($mergedIntervalsComposition ne "") {
    $mergedIntervalsFile = "logFileName = \"$mergedIntervalsComposition\"";
  }

  # Write the script file
  open(OUT, ">${scriptFile}");
  print OUT << "Header";
source("${overlappingScript}")
source("${syntenyK2Script}")
source("${breakpointsScript}")
source("${utilScript}")
source("${intervalsScript}")

print("Reading orthology table...")
orthologyTable = readOrthologyTableFile(orthologyTableFile="${tableFile}", hasHeader=F)

if (!hasRepeatedGenes(orthologyTable=orthologyTable)) {

Header

  if ($type eq "B" && $chromosomeFileGR ne "") {
    print OUT "  print(\"Reading chromosome boundaries table (genome GO)...\")\n";
    print OUT "  chromosomeBoundariesGR = readChromosomeBoundariesFile(chromosomeBoundariesFile=\"${chromosomeFileGR}\", hasHeader=F)\n";
  } else {
    print OUT "  chromosomeBoundariesGR = NULL\n";
  }

  if ($chromosomeFileGO ne "") {
    print OUT "  print(\"Reading chromosome boundaries table (genome GR)...\")\n";
    print OUT "  chromosomeBoundariesGO = readChromosomeBoundariesFile(chromosomeBoundariesFile=\"${chromosomeFileGO}\", hasHeader=F)\n";
  } else {
    print OUT "  chromosomeBoundariesGO = NULL\n";
  }

  if ($centromereFile ne "") {
    print OUT "  print(\"Reading centromere boundaries table...\")\n";
    print OUT "  centromereBoundaries = readCentromereBoundariesFile(\"${centromereFile}\", hasHeader=F)\n\n";
  }

  if ($type eq "G") {

    print OUT "  print(\"Filtering overlapping genes...\")\n";
    print OUT "  withoutOverlapping = filterOverlappingGenes(orthologyTable=orthologyTable, ${mergedIntervalsFile})\n\n";

    if ($nonOverlappingIntervalsFile ne "") {
      print OUT "  write.table(withoutOverlapping\$nonOverlappingGenesAndMergedIntervals, \"${nonOverlappingIntervalsFile}\", sep=\"\t\", quote=F, col.names=T, row.names=F)\n\n";
    }

    print OUT "  print(\"Removing type I and type II conflicts...\")\n";
    print OUT "  filteredOrthologyTable = filterOrthologyTableK2(orthologyTable=withoutOverlapping\$nonOverlappingGenesAndMergedIntervals)\n\n";

    if ($syntenyBlocksFile ne "") {
      print OUT "  write.table(filteredOrthologyTable, \"${syntenyBlocksFile}\", sep=\"\t\", quote=F, col.names=T, row.names=F)\n\n";
    }

  }

  print OUT "  print(\"Identifying breakpoints...\")\n";

  if ($type eq "G") {
    # Genes
    print OUT "  breakpoints = identifyBreakpoints1(orthologyTable=filteredOrthologyTable, limitR=${limitR}, limitAB=${limitAB}, chromosomeBoundariesGO=chromosomeBoundariesGO, extended=${extended}, extendedAB=${extendedAB}, extendBeforeVerifyLength=${extendBeforeVerifyLength}, minimumSequenceSize=${minimumSequenceSize})\n\n";
  } else {
    # Blocks
    print OUT "  breakpoints = identifyBreakpoints2(orthologyTable=orthologyTable, limitR=${limitR}, limitAB=${limitAB}, chromosomeBoundariesGR=chromosomeBoundariesGR, chromosomeBoundariesGO=chromosomeBoundariesGO, extended=${extended}, extendedAB=${extendedAB}, extendBeforeVerifyLength=${extendBeforeVerifyLength}, minimumSequenceSize=${minimumSequenceSize})\n\n";
  }

  print OUT "  print(paste(\"Number of breakpoints: \", nrow(breakpoints), sep=\"\"))\n";

  if ($centromereFile ne "") {
    print OUT "  breakpoints = removeCentromereBreakpoints(breakpoints=breakpoints, centromereBoundaries=centromereBoundaries)\n";
    print OUT "  print(paste(\"Number of breakpoints without centromere: \", nrow(breakpoints)))\n";
  }

  print OUT << "Bottom";

  print("Writing breakpoints file...")
  write.table(breakpoints, "${outputFile}", sep="\\t", quote=F, col.names=F, row.names=F)

} else {
  print("Warning: The orthology table has duplicated genes.")
}
Bottom

  close(OUT);

  # Execute R with the script
  my $logR = "$outputFile.log.Rout";
  system("$R_CMD $scriptFile $logR");
  deleteFile("$scriptFile");

}
