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
# Process the parameters
my %PARAMETERS;
verifyConfiguration();

# Mask the sequence
runRepeatMasker();

exit(0);

###############################################################################
## AUXILIARY SUB-ROUTINES #####################################################
###############################################################################

###############################################################################
# This function verifies the configuration
sub verifyConfiguration {

  my $nParameters = @ARGV;

  # Verify the number of parameters
  if ($nParameters % 2 != 0) {
    printUsage("Invalid number of parameters");
  }

  # Read the parameters
  for (my $i = 0; $i < $nParameters; $i +=2) {

    my $parameter = $ARGV[$i];
    my $value = $ARGV[$i + 1];

    $parameter =~ s/^\s+//g;
    $parameter =~ s/\s+$//g;
    $value =~ s/^\s+//g;
    $value =~ s/\s+$//g;

    if ($parameter !~ /^\-[iolsrS]$/) {
      printUsage("Unknow parameter: $parameter")
    }

    if ($value eq "") {
      printUsage("Empty string for the parameter $parameter");
    }

    $PARAMETERS{$parameter} = $value;
  }

  # Verify parameters
  # FASTA
  if (! defined $PARAMETERS{"-i"}) {
    printUsage("Missing parameter -i");
  } else {
    $PARAMETERS{"-i"} = fileExists($PARAMETERS{"-i"});
  }

  # Output file
  if (! defined $PARAMETERS{"-o"}) {
    printUsage("Missing parameter -o");
  }
  $PARAMETERS{"-o"} = getAbsolutePath($PARAMETERS{"-o"});

  # Parameter l
  if (! defined $PARAMETERS{"-l"}) {
    printUsage("Missing parameter -l");
  } else {
    if ($PARAMETERS{"-l"} ne "1" && $PARAMETERS{"-l"} ne "0") {
      printUsage("Invalid value for the parameter -l");
    }
  }

  # Parameter s
  if (! defined $PARAMETERS{"-s"}) {
    printUsage("Missing parameter -s");
  } else {
    if ($PARAMETERS{"-s"} ne "1" && $PARAMETERS{"-s"} ne "0") {
      printUsage("Invalid value for the parameter -s");
    }
  }

  # RepeatMasker path
  if (! defined $PARAMETERS{"-r"}) {
    $PARAMETERS{"-r"} = "RepeatMasker";
  } else {
    my $rm1 = getAbsolutePath($PARAMETERS{"-r"});
    if (! -e $rm1) {
      printUsage("Could not find $rm1");
    } else {
      if (-d $rm1) {
	my $rm2 = "$rm1/RepeatMasker";
	if (! -e $rm2) {
	  printUsage("$rm1 is a directory. Could not find $rm2");
	} else {
	  $PARAMETERS{"-r"} = $rm2;
	}
      }
    }
    my $rm = getAbsolutePath($PARAMETERS{"-r"});
    if (! -x $rm) {
      printUsage("$rm cannot be executed");
    }
  }

  # Parameter S
  if (! defined $PARAMETERS{"-S"}) {
    $PARAMETERS{"-S"} = "";
  } else {
    my $aux = $PARAMETERS{"-S"};
    $aux =~ s/\"//g;
    $aux =~ s/^\s+//g;
    $aux =~ s/\s+$//g;
    $aux =~ s/\s/ /g;
    if ($aux =~ / /) {
      $aux = "\"$aux\"";
    }
    $PARAMETERS{"-S"} = $aux;
  }

}

###############################################################################
# This function prints the script usage mode
sub printUsage {

  my ($message) = @_;

  if (defined $message) {
    print STDERR "\nERROR: $message\n";
  }

  print << "End_Print_Usage";

Usage:
  perl maskBreakpointSequence.pl -i fasta -o output -l [1/0] -s [1/0] -r path -S species

-------------------------------------------------------------------------------
Mandatory parameters:

  -i fasta

     "fasta" is the FASTA file that has the sequence that will be masked.

  -o output

     "output" is the name of the file that will receive the masked sequence.

  -l [1/0]

     Low complexity masking - If set with the value 1, the script will mask
     just the low complexity sequences (option -noint of the RepeatMasker).
     If set with the value 0, the script will perform full repeat masking.

  -s [1/0]

     Save RepeatMasker output files - If set with the value 1, the script
     will save the files .out and .tbl, which are generated by the
     RepeatMasker. These files will receive the names "output".out and
     "output".tbl (see parameter -o).

-------------------------------------------------------------------------------
Optional parameters:

  -r rpath

     "rpath" is the path for the program RepeatMasker. It this parameter is
     not present, the script will assume that the path of the program is listed
     in the environment variable PATH (no additional verification will be
     performed). Example: -r "/usr/local/RepeatMasker/RepeatMasker"

  -S species
 
    "species" is the value of the parameter -species of the program RepeatMasker.
    For example, -S "ciona savignyi" or -S human.

End_Print_Usage

  if (defined $message) {
    print "\nERROR: $message\n\n";
  }

  exit(1);
}

###############################################################################
# This function runs the RepeatMasker for a given FASTA file
sub runRepeatMasker {

  my $repeatmasker = $PARAMETERS{"-r"};
  my $fasta = $PARAMETERS{"-i"};
  my $low = $PARAMETERS{"-l"};
  my $save = $PARAMETERS{"-s"};
  my $output = $PARAMETERS{"-o"};
  my $species = $PARAMETERS{"-S"};

  # Create the workdir
  my $workDir = "/tmp/".rand().time();
  `mkdir $workDir`;

  # Create a soft link for the fasta file
  my $softLink = "$workDir/".rand().time();
  `ln -s $fasta $softLink`;

  # Create the command line
  my $command = "$repeatmasker -qq ";
  if ($low eq "1") {
    $command .= " -noint ";
  }
  if ($species ne "") {
    $command .= " -species $species ";
  }
  $command .= " -dir $workDir $softLink > /dev/null";
  `$command`;

  # Put the result in the right place
  if (-e "${softLink}.masked") {
    `cp -f ${softLink}.masked $output`;
  } else {
    `cp -f $fasta $output`;
  }

  if ($save eq "1") {
    if (-e "${softLink}.out") {
      `cp -f ${softLink}.out ${output}.out`;
    }
    if (-e "${softLink}.tbl") {
      `cp -f ${softLink}.tbl ${output}.tbl`;
    }
  }

  # Clean the mess
  `rm -fr $workDir`;
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
