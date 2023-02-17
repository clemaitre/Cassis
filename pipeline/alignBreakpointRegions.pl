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

# Identify the breakpoint files according with the configuration
my %BREAKPOINT_FILES;
getBreakpointRegionsFiles();

# Align the sequences
alignSequences();

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
  perl alignBreakpointRegions.pl <fastadir> <aligndir> <lastzparameters> <matrix> <lastz> [tmpdir]

-------------------------------------------------------------------------------
Mandatory parameters:

  <fastadir> Path for the directory where the FASTA file of the sequences
             SR, SA and, SB can be found.

  <aligndir> Path for the directory where the script will write the results
             of the alignment of the sequences SR against SA and SR against
             SB.

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

-------------------------------------------------------------------------------
File formats:

  + FASTA files

  The FASTA files of the sequences SR, SA and, SB must have just one sequence
  per file. The name of the files must have the following pattern:
  Breakpoint_[N]_S[X].fasta, where [N] is the number of the breakpoint and 
  [X] is one of the letters R, A or, B.

-------------------------------------------------------------------------------
Script output:

The script will align the sequences SR against SA and SR againt SB and write
the results in the directory <aligndir>. The name of the result files have the 
following pattern: Breakpoint_[N]_S[X].lastz, where [N] is the number of the 
breakpoint and [X] is one of the letters A or B.

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
  
  if ($nParameters < 5 || $nParameters > 6) {
    printUsage("Invalid number of parameters.");
  }
  
  my $fastadir    = directoryExists(trim($parameters[0]), \&printUsage);
  my $aligndir    = directoryExists(trim($parameters[1]), \&printUsage);
  
  my $lastzParameters = trim($parameters[2]);

  my $lastzMatrix = trim($parameters[3]);
  if ($lastzMatrix !~ /^default$/i) {
    $lastzMatrix = fileExists($lastzMatrix, \&printUsage);
  } else {
    $lastzMatrix =~ tr/A-Z/a-z/;
  }

  my $lastzPath = fileExists(trim($parameters[4]), \&printUsage);
  
  my $tmpDir   = "/tmp";
  if ($nParameters == 6) {
    $tmpDir = directoryExists(trim($parameters[5]), \&printUsage);
  }
  
  $PARAMETERS{"FASTADIR"}         = $fastadir;
  $PARAMETERS{"ALIGNDIR"}         = $aligndir;
  $PARAMETERS{"LASTZPARAMETERS"}  = $lastzParameters;
  $PARAMETERS{"LASTZMATRIX"}      = $lastzMatrix;
  $PARAMETERS{"LASTZPATH"}        = $lastzPath;
  $PARAMETERS{"TMPDIR"}           = $tmpDir;
  
}

###############################################################################
# This function verifies with the core script is ok
sub verifyEnvironment {

  my $corescript = "${Bin}/../core/alignSequences.pl";
  if (! -e $corescript) {
    printUsage("Missing script ${corescript}");
  }
  
  $ENVIRONMENT{"CORESCRIPT"} = $corescript;
  
}

###############################################################################
# This function reads the source directory to identify the fasta files
sub getBreakpointRegionsFiles {

  my $fastaDir = $PARAMETERS{"FASTADIR"};

  my @files = `ls -1 $fastaDir/`;
  chomp(@files);

  my %discarded;

  foreach my $file (@files) {
    
    if ($file =~ /^Breakpoint_(\d+)_S[RAB]\.fasta$/) {

      my $number = $1;
      my $key    = sprintf("%04d", $number);
      
      if (! defined $BREAKPOINT_FILES{$key} && ! defined $discarded{$key}) {
	
	my $sr = "$fastaDir/Breakpoint_${number}_SR.fasta";
	my $sa = "$fastaDir/Breakpoint_${number}_SA.fasta";
	my $sb = "$fastaDir/Breakpoint_${number}_SB.fasta";
	my $ok = 1;
	
	if (! -e $sr) {
	  print "Warning: Missing file $sr [discarding $file]\n";
	  $ok = 0;
	}
	if ($ok == 1 && ! -e $sa) {
	  print "Warning: Missing file $sa [discarding $file]\n";
	  $ok = 0;
	}
	if ($ok == 1 && ! -e $sb) {
	  print "Warning: Missing file $sb [discarding $file]\n";
	  $ok = 0;
	}
	
	if ($ok == 1) {
	  $BREAKPOINT_FILES{$key} = "$sr;$sa;$sb";
	} else {
	  $discarded{$key} = 1;
	}

      }
      
    }
    
  }
  
}

###############################################################################
# This function runs the program lastz (SR vs SA) and (SR vs SB)
sub alignSequences {

  # Run the LASTZ (SR against SA) and (SR against SB)
  foreach my $key (sort keys %BREAKPOINT_FILES) {
    my ($sr, $sa, $sb) = split(";", $BREAKPOINT_FILES{$key});
    print "LASTZ $key - SR vs SA\n";
    runLASTZ($sr, $sa, $key, 1);
    print "LASTZ $key - SR vs SB\n";
    runLASTZ($sr, $sb, $key, 0);
  }

}

###############################################################################
# This function calls the script that runs the program LASTZ
sub runLASTZ {

  my ($seq1, $seq2, $number, $sa) = @_;

  my $alignDir   = $PARAMETERS{"ALIGNDIR"};
  my $script     = $ENVIRONMENT{"CORESCRIPT"};
  my $parameters = $PARAMETERS{"LASTZPARAMETERS"};
  my $matrix     = $PARAMETERS{"LASTZMATRIX"};
  my $lastz      = $PARAMETERS{"LASTZPATH"};
  my $tmpDir     = $PARAMETERS{"TMPDIR"};

  my $lastzOutput;
  if ($sa == 1) {
    $lastzOutput = "${alignDir}/Breakpoint_${number}_SA.lastz";
  } else {
    $lastzOutput = "${alignDir}/Breakpoint_${number}_SB.lastz";
  }

  my $cmd = "perl ${script} ${seq1} ${seq2} ${lastzOutput} \"${parameters}\" ${matrix} ${lastz} ${tmpDir}";
  system($cmd);
}
