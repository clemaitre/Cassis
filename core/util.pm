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
# This function trims a string
sub trim {
  my ($string) = @_;
  if (defined $string) {
    $string =~ s/^\s+//g;
    $string =~ s/\s+$//g;
    return $string;
  }
  return "";
}

###############################################################################
# This function returns the absolute path of a file
sub getAbsolutePath {
  my ($relativePath) = @_;
  return File::Spec->rel2abs($relativePath);
}

###############################################################################
# This function verifies if the file exists.
# If the file exists, the fuction returns the absolute path for it, 
# otherwise, it returns an empty string.
# If the file does not exists and the parameter error function is defined,
# the function calls this function with a error message as parameter.
sub fileExists {
  my ($relativePath, $errorFunction) = @_;
  my $file = getAbsolutePath($relativePath);
  if (! -e $file) {
    if (defined $errorFunction) {
      &$errorFunction("Could not find the file ${file}");
    }
    return "";
  }
  if (! -f $file) {
    if (defined $errorFunction) {
      &$errorFunction("${file} is not a file.");
    }
    return "";
  }
  return $file;
}

###############################################################################
# This function verifies if the directory exists.
# If the directory exists, the fuction returns the absolute path for it, 
# otherwise, it returns an empty string.
# If the directory does not exists and the parameter error function is defined,
# the function calls this function with a error message as parameter.
sub directoryExists {
  my ($relativePath, $errorFunction) = @_;
  my $directory = getAbsolutePath($relativePath);
  if (! -e $directory) {
    if (defined $errorFunction) {
      &$errorFunction("Could not find the directory ${directory}");
    }
    return "";
  }
  if (! -d $directory) {
    if (defined $errorFunction) {
      &$errorFunction("${directory} is not a directory.");
    }
    return "";
  }
  return $directory;
}

###############################################################################
# This function deletes a file
sub deleteFile {
  my ($file) = @_;
  if (defined $file) {
    $file = getAbsolutePath($file);
    if (-e $file && -f $file) {
      `rm -f ${file}`;
    }
  }
}

###############################################################################
# This function deletes all files and subdirectories of a given directory
sub deleteDirectoryContent {
  my ($directory) = @_;
  if (defined $directory) {
    $directory = getAbsolutePath($directory);
    if (-e $directory && -d $directory) {
      `rm -fr ${directory}/*`;
    }
  }
}

###############################################################################
# This function creates a directory if it does not exist
sub createDirectory {
  my ($directory) = @_;
  if (defined $directory) {
    $directory = getAbsolutePath($directory);
    if (! -e $directory) {
      `mkdir ${directory}`;
    }
  }
}

1;
