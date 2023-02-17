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
# This file contains the definition of the function dotplotBreakpointRegion
#
# Author: Claire Lemaitre
# Reviewed by: Christian Baudet
# #############################################################################

# #############################################################################
# FUNCTION dotplotBreakpointRegion(breakpointId, alignmentSRSA, alignmentSRSB,
#                                  breakpointTable, maskedIntervals,
#                                  geneList)
#
# About the function:
#   This function plots a dotplot of a breakpoint region.
#
# Parameters:
#
#  + breakpointId = Number of the breakpoint that should be plotted
#
#  + alignmentSRSA = Dataframe that has the result of the alignment of the
#                    sequences SR and SA with lastz.
#
#  + alignmentSRSB = Dataframe that has the result of the alignment of the
#                    sequences SR and SB with lastz.
#
#  + breakpointTable = Dataframe that has the information about the
#                      breakpoints. It must have information about the
#                      breakpoint "breakpointId".
#
#  + maskedIntervals = Dataframe that has two columns: begin and end positions
#                      of the masked regions on the sequence SR.
#
#
#  + geneList = Dataframe that holds the list of genes of the first species
#               (from where the sequence SR was extracted)
#
# Dataframe formats:
#
#  + alignmentSRSA and alignmentSRSB
#    These two dataframes must have the following columns
#    - score = Score of the alignment
#    - b1    = Begin position on SR
#    - b2    = Begin position on SA (or SB)
#    - e1    = End position on SR
#    - e2    = End position on SA  (or SB)
#
#  + breakpointTable
#    This dataframe must have the following columns
#    - id        = Breakpoint ID
#    - type      = Type of the breakpoint: inter or intra
#    - sRgeneA   = Name of the gene A on the sequence sR (first genome)
#    - sRgeneB   = Name of the gene B on the sequence sR (first genome)
#    - sRchr     = Chromosome of the genes A and B (first genome)
#    - sRstrandA = Strand of the gene A (first genome)
#    - sRstrandB = Strand of the gene B (first genome)
#    - sRinf     = Inferior boundary of the sequence sR
#    - sRsup     = Superior boundary of the sequence sR
#    - sAgene    = Name of the gene A on the sequence sA (second genome)
#    - sBgene    = Name of the gene B on the sequence sB (second genome)
#    - sAchr     = Chromosome of the gene A (second genome)
#    - sBchr     = Chromosome of the gene B (second genome)
#    - sAstrand  = Strand of the gene A (second genome)
#    - sBstrand  = Strand of the gene B (second genome)
#    - sAinf     = Inferior boundary of the sequence sA
#    - sAsup     = Superior boundary of the sequence sA
#    - sBinf     = Inferior boundary of the sequence sB
#    - sBsup     = Superior boundary of the sequence sB
#    - bkpBegin  = Relative position of the breakpoint begin (related to sRinf)
#    - bkpEnd    = Relative position of the breakpoint end (related to sRinf)
#
#  + maskedIntervals
#    This dataframe must have the following columns:
#    - begin = Begin position of a masked interval on the sequence R
#    - end   = End position of a masked interval on the sequence R
#
#  + geneList
#    This dataframe must have the following columns:
#    - gene   = Name of the gene
#    - chr    = Chromosome where the gene is located
#    - inf    = Begin position of the gene
#    - sup    = End position of the gene
#    - strand = Strand of the gene
#
# PLOT:
#  The script will produce a dotplot of the alignments of the sequences SR
#  with SA and SR with SB. The fragments of the alignments of the sequences
#  SR with SA will receive the red color and the alignments of the sequences
#  SR with SB will receive the green color. The horizontal pink and turquoise
#  lines indicate, respectively, the end of the sequences SA and SB. If the
#  dataframe maskedIntervals is not NULL, the plot will have vertical light
#  yellow bars to indicate the positions of the sequence SR which are masked.
#  If the dataframe geneList is not NULL, the light blue bars indicate the
#  positions of the sequence SR which are inside of a gene. Finally, the
#  vertical dark blue lines indicate the begin and end positions of the
#  breakpoint region.
#
###############################################################################
dotplotBreakpointRegion = function(breakpointId, alignmentSRSA, alignmentSRSB,
                                   breakpointTable, maskedIntervals=NULL,
                                   geneList=NULL) {

  # Get the line that has the breakpoint which is identified
  # by "breakpointId"
  breakpointRegion = breakpointTable[breakpointTable$id==breakpointId,]

  # Get the begin and end position of the sequence SR
  beginR = breakpointRegion[1, 8]
  endR   = breakpointRegion[1, 9]

  # Get the begin and end position of the sequence SA
  beginA = breakpointRegion[1, 16]
  endA   = breakpointRegion[1, 17]

  # Get the begin and end position of the sequence SB
  beginB = breakpointRegion[1, 18]
  endB   = breakpointRegion[1, 19]

  # Calculate the length of the sequences SR, SA and, SB
  lengthR = endR - beginR + 1
  lengthA = endA - beginA + 1
  lengthB = endB - beginB + 1

  # Get the strand of the gene A on the sequences SR and SA
  strandGeneAonSR = breakpointRegion[1, 6]
  strandGeneAonSA = breakpointRegion[1, 14]

  # Get the strand of the gene B on the sequences SR and SB
  strandGeneBonSR = breakpointRegion[1, 7]
  strandGeneBonSB = breakpointRegion[1, 15]

  # Define the signal of the strand of the sequences SA and SB
  strandA = ifelse(strandGeneAonSR==strandGeneAonSA, "+", "-")
  strandB = ifelse(strandGeneBonSR==strandGeneBonSB, "+", "-")

  # Get the relative position of the begin and end of the breakpoint
  breakpointBegin = breakpointRegion[1, 20]
  breakpointEnd   = breakpointRegion[1, 21]

  # Get the list of genes that are inside of the interval defined
  # by the sequence SR
  filteredGeneList = NULL
  if (!is.null(geneList)) {
    chromossomeSR   = as.vector(breakpointRegion[1, 5])
    beginR          = as.vector(beginR)
    endR            = as.vector(endR)
    filteredGeneList = geneList[geneList$chr==chromossomeSR & geneList$inf>=beginR & geneList$sup<=endR, ]
    filteredGeneList$inf = filteredGeneList$inf - beginR
    filteredGeneList$sup = filteredGeneList$sup - beginR
  }

  # Define the title, xlabel and ylabel of the plot
  title = paste("Breakpoint ", breakpointId, ": SA", strandA, " SB", strandB)
  labelX = "SR"
  labelY = "SA and SB"

  doubleDotplot(alignmentSRSA=alignmentSRSA,
                alignmentSRSB=alignmentSRSB,
                lengthR=lengthR,
                lengthA=lengthA,
                lengthB=lengthB,
                breakpointBegin=breakpointBegin,
                breakpointEnd=breakpointEnd,
                title=title,
                labelX=labelX,
                labelY=labelY,
                maskedIntervals=maskedIntervals,
                filteredGeneList)
}




# #############################################################################
# FUNCTION doubleDotplot(alignmentSRSA, alignmentSRSB, lengthR, lengthA,
#                        lengthB, breakpointBegin, breakpointEnd, title,
#                        labelX, labelY, maskedIntervals, geneList)
#
# About the function:
#   This function plots a dotplot of a breakpoint region.
#
# Parameters:
#
#  + alignmentSRSA = Dataframe that has the result of the alignment of the
#                    sequences SR and SA with lastz.
#
#  + alignmentSRSB = Dataframe that has the result of the alignment of the
#                    sequences SR and SB with lastz.
#
#  + lengthR = Length of the sequence SR.
#
#  + lengthA = Length of the sequence SA.
#
#  + lengthB = Length of the sequence SB.
#
#  + breakpointBegin = Begin position of the breakpoint.
#
#  + breakpointEnd = End position of the breakpoint.
#
#  + title = Title of the plot.
#
#  + labelX = Label of X axis.
#
#  + labelY = Label of Y axis.
#
#  + maskedIntervals = Dataframe that has two columns: begin and end positions
#                      of the masked regions on the sequence SR.
#
#  + geneList = List of genes that are inside of the sequence SR.
#
# Dataframe formats:
#
#  + alignmentSRSA and alignmentSRSB
#    These two dataframes must have the following columns
#    - score = Score of the alignment
#    - b1    = Begin position on SR
#    - b2    = Begin position on SA (or SB)
#    - e1    = End position on SR
#    - e2    = End position on SA  (or SB)
#
#  + maskedIntervals
#    This dataframe must have the following columns:
#    - begin = Begin position of a masked interval on the sequence R
#    - end   = End position of a masked interval on the sequence R
#
#  + geneList
#    This dataframe must have the following columns:
#    - gene   = Name of the gene
#    - chr    = Chromosome where the gene is located
#    - inf    = Begin position of the gene
#    - sup    = End position of the gene
#    - strand = Strand of the gene
#
# PLOT:
#  The script will produce a dotplot of the alignments of the sequences SR
#  with SA and SR with SB. The fragments of the alignments of the sequences
#  SR with SA will receive the red color and the alignments of the sequences
#  SR with SB will receive the green color. The horizontal pink and turquoise
#  lines indicate, respectively, the end of the sequences SA and SB. If the
#  dataframe maskedIntervals is not NULL, the plot will have vertical light
#  yellow bars to indicate the positions of the sequence SR which are masked.
#  If the dataframe geneList is not NULL, the light blue bars indicate the
#  positions of the sequence SR which are inside of a gene. Finally, the
#  vertical dark blue lines indicate the begin and end positions of the
#  breakpoint region.
###############################################################################
doubleDotplot = function(alignmentSRSA, alignmentSRSB, lengthR, lengthA, lengthB,
                         breakpointBegin, breakpointEnd, title, labelX, labelY,
                         maskedIntervals=NULL, geneList=NULL) {

  # Define the maxY of the plot (maximum between lenghtA and lenghtB)
  maxY = max(lengthA, lengthB)

  # Create the plot area
  plot( c(0,lengthR), c(0,maxY), type="n", main=title, xlab=labelX, ylab=labelY)

  # Light blue to identify the positions where we have genes
  if (!is.null(geneList)) {
    if (nrow(geneList) > 0) {
      rect(geneList$inf, 1, geneList$sup, maxY, col = "lightsteelblue1", border=NA)
    }
  }

  # Light yellow lines to identify positions that were masked on the sequence SR
  if (!is.null(maskedIntervals)) {
    maskedIntervalsOnSR = unionIntervals(maskedIntervals)
    if (nrow(maskedIntervalsOnSR) > 0) {
      rect(maskedIntervalsOnSR$begin, 1, maskedIntervalsOnSR$end, maxY,
           col = "lightyellow2", border=NA)
    }
  }

  # Plot the fragments of the alignment between SR and SA (red segments)
  if (nrow(alignmentSRSA) > 0) {
    segments(alignmentSRSA$b1, alignmentSRSA$b2, alignmentSRSA$e1, alignmentSRSA$e2, col=2)
  }
  # Horizontal line to identify the end of the sequence SA (pink line)
  abline(h=lengthA, col="pink3")

  # Plot the fragments of the alignment between SR and SB (green segments)
  if (nrow(alignmentSRSB) > 0) {
    segments(alignmentSRSB$b1, alignmentSRSB$b2, alignmentSRSB$e1, alignmentSRSB$e2, col=3)
  }
  # Horizontal line to identify the end of the sequence SB (turquoise line)
  abline(h=lengthB, col="turquoise")

  # Vertical lines to identify the begin and the end of the breakpoint
  # (dark blue lines)
  abline(v=breakpointBegin, col=4)
  abline(v=breakpointEnd, col=4)
}

