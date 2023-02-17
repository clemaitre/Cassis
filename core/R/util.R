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

# #############################################################################
# This file contains the definition of a set of auxiliary functions
#
# Author: Christian Baudet
# #############################################################################

# #############################################################################
# FUNCTION readChromosomeBoundariesFile(chromosomeBoundariesFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that is specified by the parameter
#     chromosomeBoundariesFile. This file must have the following columns:
#     chr  - Name of the chromosome
#     inf  - Leftmost extremity of the chromosome
#     sup  - Rightmost extremity of the chromosome
#     file - Name of the chromosome FASTA file
#
# 2 - The function requires that the file must have the columns in the
#     same order as described above and that the chromosomes are related to
#     just one species.
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them to "chr", "inf", "sup" and, "file".
#
# 5 - The function reads the file and returns a dataframe with the columns
#     "chr", "inf", "sup" and, "file".
#
# #############################################################################
readChromosomeBoundariesFile = function(chromosomeBoundariesFile, hasHeader=F) {
  chromosomeBoundaries = read.table(chromosomeBoundariesFile, header=hasHeader)
  names(chromosomeBoundaries) = c("chr", "inf", "sup", "file")
  return (chromosomeBoundaries)
}

# #############################################################################
# FUNCTION readOrthologyTableFile(orthologyTableFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that is specified by the parameter
#     orthologyTableFile. This file must have the following columns:
#     g1      - Name of the gene in the first genome
#     c1      - Chromosome where the gene g1 is located
#     inf1    - Start position of the gene g1 on the chromosome c1
#     sup1    - End position of the gene g1 on the chromosome c1
#     strand1 - Strand of the gene g1
#     g2      - Name of the gene in the second genome
#     c2      - Chromosome where the gene g2 is located
#     inf2    - Start position of the gene g2 on the chromosome c2
#     sup2    - End position of the gene g2 on the chromosome c2
#     strand2 - Strand of the gene g2
#
# 2 - The function requires that the file must have the columns in the
#     same order as described above.
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them with the names which are listed above.
#
# 5 - The function reads the file and returns a dataframe with the columns
#     which are listed above.
#
# #############################################################################
readOrthologyTableFile = function(orthologyTableFile, hasHeader=F) {
  orthologyTable = read.table(orthologyTableFile, header=hasHeader)
  names(orthologyTable) = c("g1", "c1", "inf1", "sup1", "strand1", "g2", "c2", "inf2", "sup2", "strand2")
  return (orthologyTable)
}

# #############################################################################
# FUNCTION readCentromereBoundariesFile(centromereBoundariesFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that is specified by the parameter
#     centromereBoundariesFile. This file must have the following columns:
#     chr - Name of the chromosome
#     inf - Leftmost extremity of the centromere of the chromosome
#     sup - Rightmost extremity of the centromere of the chromosome
#
# 2 - The function requires that the file must have the columns in the
#     same order as described above and that the chromosomes are related to
#     just one species.
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them to "chr", "inf" and, "sup".
#
# 5 - The function reads the file and returns a dataframe with the columns
#     "chr", "inf" and, "sup".
#
# #############################################################################
readCentromereBoundariesFile = function(centromereBoundariesFile, hasHeader=F) {
  centromereBoundaries = read.table(centromereBoundariesFile, header=hasHeader)
  names(centromereBoundaries) = c("chr", "inf", "sup")
  return (centromereBoundaries)
}

# #############################################################################
# FUNCTION readGeneListFile(geneListFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that is specified by the parameter
#     geneListFile. This file must have the following columns:
#     gene   - Name of the gene
#     chr    - Chromosome where the gene is located
#     inf    - Start position of the gene the chromosome
#     sup    - End position of the gene on the chromosome
#     strand - Strand of the gene
#
# 2 - The function requires that the file must have the columns in the
#     same order as described above.
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them with the names which are listed above.
#
# 5 - The function reads the file and returns a dataframe with the columns
#     which are listed above.
#
# #############################################################################
readGeneListFile = function(geneListFile, hasHeader=F) {
  geneList = read.table(geneListFile, header=hasHeader)
  names(geneList) = c("gene", "chr", "inf", "sup", "strand")
  return (geneList)
}

# #############################################################################
# FUNCTION readBreakpointsFile(breakpointsFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that is specified by the parameter
#     chromosomeBoundariesFile. This file must have the following columns:
#     id        - Breakpoint ID
#     type      - Type of the breakpoint: inter or intra
#     sRgeneA   - Name of the gene A on the sequence sR (first genome)
#     sRgeneB   - Name of the gene B on the sequence sR (first genome)
#     sRchr     - Chromosome of the genes A and B (first genome)
#     sRstrandA - Strand of the gene A (first genome)
#     sRstrandB - Strand of the gene B (first genome)
#     sRinf     - Inferior boundary of the sequence sR
#     sRsup     - Superior boundary of the sequence sR
#     sAgene    - Name of the gene A on the sequence sA (second genome)
#     sBgene    - Name of the gene B on the sequence sB (second genome)
#     sAchr     - Chromosome of the gene A (second genome)
#     sBchr     - Chromosome of the gene B (second genome)
#     sAstrand  - Strand of the gene A (second genome)
#     sBstrand  - Strand of the gene B (second genome)
#     sAinf     - Inferior boundary of the sequence sA
#     sAsup     - Superior boundary of the sequence sA
#     sBinf     - Inferior boundary of the sequence sB
#     sBsup     - Superior boundary of the sequence sB
#     bkpBegin  - Relative position of the breakpoint begin (related to sRinf)
#     bkpEnd    - Relative position of the breakpoint end (related to sRinf)
#     status    - Breakpoint status
#
# 2 - The function requires that the file must have the columns in the
#     same order as described above.
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them with the names which are listed above.
#
# 5 - The function reads the file and returns a dataframe which have columns
#     that have the names which are listed above.
#
# #############################################################################
readBreakpointsFile = function(breakpointsFile, hasHeader=F) {
  breakpoints = read.table(breakpointsFile, header=hasHeader)
  names(breakpoints) = c("id", "type", "sRgeneA", "sRgeneB", "sRchr",
                         "sRstrandA", "sRstrandB", "sRinf", "sRsup",
                         "sAgene", "sBgene", "sAchr", "sBchr", "sAstrand",
                         "sBstrand", "sAinf", "sAsup", "sBinf", "sBsup",
                         "bkpBegin", "bkpEnd", "status")
  return (breakpoints)
}

# #############################################################################
# FUNCTION readLastzOutputFile(lastzOutputFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that is specified by the parameter
#     lastzOutputFile. This file must have the following columns:
#     score - Score of the alignment
#     b1    - Begin position on SR
#     b2    - Begin position on SA (or SB)
#     e1    - End position on SR
#     e2    - End position on SA  (or SB)
#
# 2 - The function requires that the file must have the columns in the
#     same order as described above.
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them with the names which are listed above.
#
# 5 - The function reads the file and returns a dataframe with the columns
#     which are listed above.
#
# #############################################################################
readLastzOutputFile = function(lastzOutputFile, hasHeader=F) {
  lastzOutput = read.table(lastzOutputFile, header=hasHeader)
  names(lastzOutput) = c("score", "b1", "b2", "e1", "e2")
  return (lastzOutput)
}

# #############################################################################
# FUNCTION hasDuplicatedGenes(orthologyTable)
#
# About this function:
# 1 - This function verifies if there are duplicated genes on the orthology
#     table received by parameter
#
# 2 - This function receives a table of ortholog genes which has the columns:
#     g1      - Name of the gene in the first genome
#     g2      - Name of the gene in the second genome
#
# 3 - The function returns TRUE if it had find duplicated genes and FALSE
#     otherwise.
#
# #############################################################################
hasRepeatedGenes = function (orthologyTable) {

  if (!is.null(orthologyTable)) {
    numberOfRows = nrow(orthologyTable)
    uniquegenes1 = length(unique(orthologyTable$g1))
    uniquegenes2 = length(unique(orthologyTable$g2))
    if (numberOfRows != uniquegenes1 | numberOfRows != uniquegenes2) {
      return (T)
    }
  }
  return (F)
}

# #############################################################################
# FUNCTION getChromosomeBegin(chromosome, chromosomeBoundaries)
#
# About this function:
#
# 1 - This function returns the begin position of the chromosome identified
#     by the given parameter.
#
# #############################################################################
getChromosomeBegin = function(chromosome, chromosomeBoundaries) {
  return (chromosomeBoundaries$inf[chromosomeBoundaries$chr == chromosome])
}

# #############################################################################
# FUNCTION getChromosomeEnd(chromosome, chromosomeBoundaries)
#
# About this function:
#
# 1 - This function returns the end position of the chromosome identified
#     by the given parameter.
#
# #############################################################################
getChromosomeEnd = function(chromosome, chromosomeBoundaries) {
  return (chromosomeBoundaries$sup[chromosomeBoundaries$chr == chromosome])
}

# #############################################################################
# FUNCTION validateBeginPosition(position, chromosome, chromosomeBoundaries)
#
# About this function:
#
# 1 - This function verifies if the given position respects the left boundary
#     of the chromosome. If the position does not respect, the method tries to
#     correct the position.
#
# #############################################################################
validateBeginPosition = function(position, chromosome, chromosomeBoundaries) {  
  if (!is.null(chromosomeBoundaries)) {
    begin = getChromosomeBegin(chromosome, chromosomeBoundaries)
    if (position < begin) {
      position = begin
    }
  } else {
    if (position < 1) {
      position = 1
    }
  }
  return (position)
}

# #############################################################################
# FUNCTION validateEndPosition(position, chromosome, chromosomeBoundaries)
#
# About this function:
#
# 1 - This function verifies if the given position respects the right boundary
#     of the chromosome. If the position does not respect, the method tries to
#     correct the position.
#
# #############################################################################
validateEndPosition = function(position, chromosome, chromosomeBoundaries) {
  if (!is.null(chromosomeBoundaries)) {
    end = getChromosomeEnd(chromosome, chromosomeBoundaries)
    if (position > end) {
      position = end
    }
  }
  return (position)
}

# #############################################################################
# FUNCTION removeCentromereBreakpoints(breakpoints, centromereBoundaries)
#
# About this function:
#
# 1 - This function returns a new data frame of breakpoints without the
#     breakpoints that have a centromere inside of them.
#
# 2 - The table breakpoints must have at least these columns:
#     id    - Breakpoint ID
#     sRchr - Chromosome of the genes A and B (first genome)
#     sRinf - Inferior boundary of the sequence sR
#     sRsup - Superior boundary of the sequence sR
#
# 3 - The table centromereBoundaries must have the columns:
#     chr - Name of the chromosome
#     inf - Leftmost extremity of the centromere of the chromosome
#     sup - Rightmost extremity of the centromere of the chromosome
#
# #############################################################################
removeCentromereBreakpoints = function(breakpoints, centromereBoundaries) {

  # Get the breakpoints that have chromosomes listed on the
  # centromereBoundaries table
  withCentromere = breakpoints[is.element(breakpoints$sRchr, centromereBoundaries$chr),]

  # Process each line of the table
  centromereInside = NULL
  for (i in 1:nrow(withCentromere)) {
    # Get the right centromere
    centromere = centromereBoundaries[as.character(centromereBoundaries$chr) == as.character(withCentromere$sRchr[i]),]
    # Verify the boundaries
    if (withCentromere$sRinf[i] < centromere$inf & withCentromere$sRsup[i] > centromere$sup) {
      # In this case, the breakpoint has a centromere inside of it
      centromereInside = c(centromereInside, withCentromere$id[i])
    }
  }

  # Remove the breakpoints which have a centromere inside of it.
  toReturn = breakpoints[!is.element(breakpoints$id, centromereInside),]
  return (toReturn)
}

# #############################################################################
# FUNCTION readRepeatMaskerOutput(repeatMaskerOutputFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that contains the formated output of the
#     program RepeatMasker (version 3.2.8 or later). To format the output of
#     the program RepeatMasker (file .out) use the script Perl called
#     formatRepeatMaskerOutput.pl that is on the directory util.
#
# 2 - The formated output is a TSV file that has the following colums:
#     SW    - Smith-Waterman score of the match
#     %sub  - % substitutions in matching region compared to the consensus
#     %gq   - % of bases opposite a gap in the query sequence (deleted bp)
#     %gr   - % of bases opposite a gap in the repeat consensus (inserted bp)
#     qname - Name of query sequence
#     qb    - Starting position of match in query sequence
#     qe    - Ending position of match in query sequence
#     qleft - No. of bases in query sequence past the ending position of match
#     compl - C = Match is with the Complement of the consensus sequence in
#             the database
#     rname - Name of the matching interspersed repeat
#     rtype - The class of the repeat
#     rb    - Starting position of match in database sequence (using
#             top-strand numbering)
#     re    - Ending position of match in database sequence
#     rleft - no. of bases in (complement of) the repeat consensus sequence
#             prior to beginning of the match
#     id    - Identification number of the match
#     inc   - An asterisk (*) in this column indicates that there is a
#             higher-scoring match whose domain partly (<80%) includes the
#             domain of this match.
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them to "sw", "sub", "gq", "gr", "qname", "qb", "qe",
#     "qleft", "compl", "rname", "rtype", "rb", "re", "rleft", "id", and "inc".
#
# 5 - The function reads the file and returns a dataframe with the columns
#     "sw", "sub", "gq", "gr", "qname", "qb", "qe", "qleft", "compl", "rname",
#     "rtype", "rb", "re", "rleft", "id", and "inc".
#
# #############################################################################
readRepeatMaskerOutput = function(repeatMaskerOutputFile, hasHeader=F) {
  repeatMaskerOutput = read.table(repeatMaskerOutputFile, header=hasHeader)
  names(repeatMaskerOutput) =
    c("sw", "sub", "gq", "gr", "qname", "qb", "qe", "qleft", "compl", "rname",
      "rtype", "rb", "re", "rleft", "id", "inc")
  return (repeatMaskerOutput)
}

# #############################################################################
# FUNCTION getMaskedIntervals(repeatMaskerOutput)
#
# About this function:
# 1 - This function get the dataframe that has the RepeatMasker output and
#     returns a dataframe that has just the begin and end positions the
#     masked intervals on the query sequence.
#
# 2 - This function expects that the input dataframe has the columns qb and qe
#     (see function readRepeatMaskerOutput)
#
# 3 - The function returns a dataframe that have the columns begin and end that
#     identify the coordinates of the masked intervals on the query sequence.
#
# #############################################################################
getMaskedIntervals = function(repeatMaskerOutput) {
  if (is.null(repeatMaskerOutput)) {
    return (NULL)
  }
  return (data.frame(begin=repeatMaskerOutput$qb, end=repeatMaskerOutput$qe))
}

# #############################################################################
# FUNCTION readMaskedIntervalsFile(maskedIntervalsFile, hasHeader)
#
# About this function:
# 1 - This function reads the file that contains the list of masked intervals
#     that exists on a DNA sequence.
#
# 2 - The file must have the following columns:
#     begin - Begin position of the masked interval
#     end   - End position of the masked interval
#
# 3 - The user must set the boolean parameter hasHeader to inform if the file
#     has or not a header row which identifies the columns (default value = F).
#
# 4 - Note that, if the file has a header, the function ignores the column
#     names and rename them to "begin" and "end".
#
# 5 - The function reads the file and returns a dataframe with the columns
#     "begin" and "end".
#
# #############################################################################
readMaskedIntervalsFile = function(maskedIntervalsFile, hasHeader=F) {
  maskedIntervals = read.table(maskedIntervalsFile, header=hasHeader)
  names(maskedIntervals) = c("begin", "end")
  return (maskedIntervals)
}
