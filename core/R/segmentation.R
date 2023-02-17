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
# This file contains the definition of the function segmentAndPlotABreak
#
# Author: Claire Lemaitre
# Reviewed by: Christian Baudet
# #############################################################################

# #############################################################################
# FUNCTION segmentAndPlotABreak(breakpointId, alignmentSRSA, alignmentSRSB,
#                               breakpointTable,  maskedIntervals,
#                               toplot)
#
# About the function:
#   This function performs the segmentation of the breakpoint region identified
#   by the parameter "breakpointId". If the parameter toplot is TRUE, the
#   function also plots the segmentation result.
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
#  + toplot = If TRUE, the function will plot the segmentation result.
#
# Function return:
#  The function will return a dataframe that has one row and the following
#  columns:
#  - id             = Breakpoint identifier.
#  - chrR           = Chromosome where the sequence R is located.
#  - srLength       = Length of the sequence R.
#  - masked         = Number of masked bases on the sequence R.
#  - hits           = Number of bases of the sequence R that have at least
#                     one hit with the sequences A or B.
#  - hitsSA         = Number of bases of the sequence R that have a hit with
#                     the sequence A
#  - hitsSB         = Number of bases of the sequence R that have a hit with
#                     the sequence B
#  - hitsSAB        = Number of bases of the sequence R that have hits with
#                     the sequence A and the sequence B
#  - oldBkpBegin    = Begin position of the original breakpoint (non refined)
#  - oldBkpEnd      = End position of the original breakpoint (non refined)
#  - oldBkpLength   = Lenght of the original breakpoint (non refined)
#  - status         = Breakpoint status (more information bellow)
#  - left           = Value of the score curve in the left extremity (begin)
#                     of the new breakpoint region.
#  - right          = Value of the score curve in the right extremity (end)
#                     of the new breakpoint region.
#  - nonMaskedLeft  = Value of the score curve in the left extremity (begin)
#                     of the new breakpoint region (discarding the masked
#                     bases).
#  - nonMaskedRight = Value of the score curve in the right extremity (end)
#                     of the new breakpoint region (discarding the masked
#                     bases).
#  - bkpBegin       = Begin position of the refined breakpoint.
#  - bkpEnd         = End position of the refined breakpoint.
#  - bkpLength      = Length of the refined breakpoint.
#  - bkpNonMasked   = Length of the refined breakpoint (discarding the masked
#                     bases).
#  - bkpHits        = Number of bases of the refined breakpoint that are
#                     covered by a lastz hit.
#  - bkpDoubleHits  = Number of bases of the refined breakpoint that are
#                     covered by a lastz hit with the sequences A and B.
#  - RSS            = RSS value of the segmentation.
#  - cbhBkpBegin    = Begin position of the refined breakpoint on the
#                     sequence R' (sequence that has just the positions
#                     that are covered by lastz hits).
#  - cbhBkpEnd      = End position of the refined breakpoint on the sequence
#                     R' (sequence that has just the positions that are 
#                     covered by lastz hits).
#  - cbhLeft        = Value of the score curve in the left extremity (begin)
#                     of the refined breakpoint region on the sequence R'.
#  - cbhRight       = Value of the score curve in the right extremity (end)
#                     of the refined breakpoint region on the sequence R'.
#
#  Breakpoint status table:
#     1 : Segmentation result passed on the statistical test 
#     0 : Segmentation result failed on the statistical test
#    -1 : Segmentation was not performed because there are no hits on
#         the alignments of the sequences SR vs SA and SR vs SB
#    -2 : Segmentation was not performed because sequence SR is smaller
#         than the allowed limit
#    -3 : Segmentation was not performed because sequence SA is smaller 
#         than the allowed limit
#    -4 : Segmentation was not performed because sequence SB is smaller 
#         than the allowed limit
#    -5 : Segmentation was not performed because sequence SA and SB are 
#         smaller than the allowed limit
#    -6 : Segmentation was not performed because sequence SR is bigger 
#         than the allowed limit
#
# #############################################################################
segmentAndPlotABreak = function(breakpointId, alignmentSRSA, alignmentSRSB,
                                breakpointTable,  maskedIntervals=NULL,
                                toplot=T) {

  # Get the line that has the breakpoint which is identified
  # by "breakpointId"
  breakpointRegion = breakpointTable[breakpointTable$id==breakpointId,]

  # Get the begin and end position of the sequence SR
  beginR = breakpointRegion[1, 8]
  endR   = breakpointRegion[1, 9]

  # Get the chromosome of the sequence SR
  chrR = breakpointRegion[1, 5]

  # Calculate the length of the sequences SR
  lengthR = endR - beginR + 1

  # Get the relative position of the begin and end of the breakpoint
  breakpointBegin = breakpointRegion[1, 20]
  breakpointEnd   = breakpointRegion[1, 21]

  # Get the breakpoint status
  status =  breakpointRegion[1, 22]


  # Read the sequence SR to get information about the masked positions
  maskedPositionsOnSR = vector()
  maskedIntervalsOnSR = data.frame()
  numberOfMaskedBases = 0
  if (status == 1 & !is.null(maskedIntervals)) {
    # Make the union of the intervals
    maskedIntervalsOnSR = unionIntervals(maskedIntervals)
    # Create the vector of masked positions
    maskedPositionsOnSR = unique(unlist(apply(maskedIntervalsOnSR, 1,
      function(x) seq(as.numeric(x[1]), as.numeric(x[2]), by = 1)),
      use.names = FALSE))
    # Calculate the number of masked bases
    nMaskedIntervals = nrow(maskedIntervalsOnSR)
    for (i in 1:nMaskedIntervals) {
      b = maskedIntervalsOnSR[i, 1]
      e = maskedIntervalsOnSR[i, 2]
      numberOfMaskedBases = numberOfMaskedBases + (e - b + 1)
    }
  }

  # Read the alignment files to access the information about the hits
  # of the alignment SR vs SA
  hitsSRSA = getHitsInformation(alignment = alignmentSRSA, lengthSeqR = lengthR)
  # of the alignment SR vs SB
  hitsSRSB = getHitsInformation(alignment = alignmentSRSB, lengthSeqR = lengthR)

  if (status == 1 & (hitsSRSA$numberOfHits > 0 || hitsSRSB$numberOfHits > 0)) {
  
    # Calculate the difference between the list of positions that have hits
    # on the alignments of SR-SA and SR-SB.
    # If diffHitsSASB[i]
    # =  1, there is a hit on the position i on the alignment SR-SA and no hit
    #       on the position i on the alignment SR-SB
    # =  0, there is a hit on the position i on the alignment SR-SA and on the
    #       alignment SR-SB
    # = -1, there is a hit on the position i on the alignment SR-SB and no hit
    #       on the position i on the alignment SR-SA
    diffHitsSASB = hitsSRSA$positions - hitsSRSB$positions

    # Clean the vector diffHitsSASB. The masked positions receive the value 0.
    diffHitsSASB[maskedPositionsOnSR] = 0

    # Create the score vector based of the diff vector and the vectors of hit
    # positions. Positions that do not show hit on the alignments of SR-SA and
    # SR-SB will receive the value 2. Otherwise, the position will receive
    # the corresponding value of the vector diffHitsSASB
    # If scoreVector[i]
    # =  2, there is no hit on the position i on the alignment SR-SA and on the
    #       alignment SR-SB
    # =  1, there is a hit on the position i on the alignment SR-SA and no hit
    #       on the position i on the alignment SR-SB
    # =  0, there is a hit on the position i on the alignment SR-SA and on the
    #       alignment SR-SB
    # = -1, there is a hit on the position i on the alignment SR-SB and no hit
    #       on the position i on the alignment SR-SA
    scoreVector = ifelse(hitsSRSA$positions == 0 & hitsSRSB$positions == 0, 2, diffHitsSASB)

    # Compute the lengths and values of runs of equal values on the scoreVector
    scoreVectorRuns = rle(scoreVector)
  
    # Create a data frame that have two columns:
    # x = position on the sequence SR
    # y = value of the position x on the vector diffHitsSASB
    scoreTable = data.frame(x = seq(1, lengthR, by = 1), y = diffHitsSASB)

    # Data frame that holds the begin and end positions of all intervals that
    # are not covered by any hit (consecutive elements that have the value 2 on
    # vector scoreVector)
    nonCoveredByHitsIntervals = getNonCoveredByHitsIntervals(scoreVectorRuns)
  
    # Create a score vector that does not have the positions that do not show
    # hits in both alignments (positions on the scoreVector that have value 2)
    scoreVectorCoveredByHits = scoreVector[scoreVector != 2]

    # Calculate the size of the sequence SR after removing the positions that
    # have no hits in both alignments
    lengthRCoveredByHits = length(scoreVectorCoveredByHits)

    # Create a data frame that have two columns:
    # x = position on the sequence SR after removing the positions that do not
    #     show hits in both alignments
    # y = value of the position x on the vector scoreVectorCoveredByHits
    dataCoveredByHits =
      data.frame(x = seq(1, lengthRCoveredByHits, by=1), y = scoreVectorCoveredByHits)

    # Create a vector that has the "cumulated length" just of the intervals that
    # are covered by hits
    coveredByHitsCumulatedIntervalLengths = getCoveredByHitsCumulatedIntervalLengths(scoreVectorRuns)

    # Perform the segmentation of the breakpoint region
    # Define the breakpoint region with a better resolution
    # (discover the position where the score curve has a plateau)
    segmentation = segmentMean(scoreVectorCoveredByHits, coveredByHitsCumulatedIntervalLengths)

    # Do the statistical test to validate the segmentation
    ok = statisticalTest(scoreVectorRuns, segmentation$RSS)
  
    
    # Get the extremities of the breakpoint on the sequence SR'
    coveredByHitsBreakpointBegin = segmentation$begin
    coveredByHitsBreakpointEnd = segmentation$end

    # Correct the extremities (map them to the sequence SR)
    newBreakpointBegin = getTruePosition(coveredByHitsBreakpointBegin, nonCoveredByHitsIntervals)
    newBreakpointEnd = getTruePosition(coveredByHitsBreakpointEnd, nonCoveredByHitsIntervals)

    # Calculate the left coefficients
    coefLeft = 0
    coefLeftNonMasked = 0
    coefLeftCoveredByHits = 0
    if (newBreakpointBegin > 0) {
      coefLeft = signif(mean(diffHitsSASB[1:newBreakpointBegin]), 3)
      coefLeftNonMasked =
        signif(sum(diffHitsSASB[1:newBreakpointBegin]) /
               (newBreakpointBegin - sum(maskedPositionsOnSR <= newBreakpointBegin)), 3)
      coefLeftCoveredByHits =
        signif(mean(scoreVectorCoveredByHits[1:coveredByHitsBreakpointBegin]), 3)
    }

    # Calculate the right coefficients
    coefRight = 0
    coefRightNonMasked = 0
    coefRightCoveredByHits = 0
    if (newBreakpointEnd < length(diffHitsSASB)) {
      coefRight = signif(mean(diffHitsSASB[newBreakpointEnd:length(diffHitsSASB)]), 3)
      coefRightNonMasked =
        signif(sum(diffHitsSASB[newBreakpointEnd:length(diffHitsSASB)]) /
               (length(diffHitsSASB) - newBreakpointEnd + 1 - sum(maskedPositionsOnSR >= newBreakpointEnd)), 3)
      coefRightCoveredByHits =
        signif(mean(scoreVectorCoveredByHits[coveredByHitsBreakpointEnd:length(scoreVectorCoveredByHits)]), 3)
    }

    if (toplot) {
      plotTitle =
        paste(paste("Breakpoint ", breakpointId, "  chr = ", chrR, sep=""),
              paste("left = ", coefLeft, "  right = ", coefRight, sep=""),
              paste("begin = ", newBreakpointBegin, "  end = ", newBreakpointEnd, sep=""),
              sep="\n")
    
      coveredByHitsPlotTitle =
        paste(paste("Breakpoint ", breakpointId, "  chr = ", chrR, sep=""),
              paste("left = ", coefLeftCoveredByHits, "  right = ", coefRightCoveredByHits, sep=""),
              paste("begin = ", coveredByHitsBreakpointBegin, "  end = ", coveredByHitsBreakpointEnd, sep=""),
              sep="\n")
    
      par(mfrow = c(1, 2))
    
      # Plot of the breakpoint region (complete)
      plotSegmentation(scoreTable,
                       hitsSRSA$intervals, hitsSRSB$intervals,
                       newBreakpointBegin, newBreakpointEnd,
                       plotTitle,
                       breakpointBegin, breakpointEnd,
                       maskedIntervalsOnSR)

      # Plot of the breakpoint region after removing the positions
      # that are not covered by hits
      plot(dataCoveredByHits$x,
           cumsum(dataCoveredByHits$y),
           type = "l",
           main = coveredByHitsPlotTitle,
           xlab = "Position on SR'",
           ylab = "Cumulated Score")
      # Draw a red and a green line to indicate the new breakpoint extremities
      abline(v = coveredByHitsBreakpointBegin, col=2, lwd=2)
      abline(v = coveredByHitsBreakpointEnd, col=3, lwd=2)
    }

    # Calculate the lenght of the breakpoint regions (complete, without masked
    # sequences and without intervals that do not have hits)
    newBreakpointLength = newBreakpointEnd - newBreakpointBegin
    coveredByHitsBreakpointLength = coveredByHitsBreakpointEnd - coveredByHitsBreakpointBegin - 1
    nonMaskedBreakpointLength = newBreakpointEnd - newBreakpointBegin -
      sum(maskedPositionsOnSR < newBreakpointEnd & maskedPositionsOnSR > newBreakpointBegin) - 1

    # Calculate the number of positions that have hit in the two alignments
    # (SR vs SA and SR vs SB)
    numberOfDoubleHitPositions = sum(scoreVectorCoveredByHits == 0)
  
    # Calculate the number of positions that have hit in the two alignments
    # (SR vs SA and SR vs SB) inside of the new breakpoint region
    numberOfDoubleHitBreakpointPositions =
      sum(scoreVectorCoveredByHits[coveredByHitsBreakpointBegin:coveredByHitsBreakpointEnd] == 0)

    # Calculate the old breakpoint length to add into the final dataframe
    oldBreakpointLength = breakpointEnd - breakpointBegin
  
    # Produce a data frame that groups all segmentation information
    segmentation =
      data.frame(id             = breakpointId,
                 chrR           = chrR,
                 srLength       = lengthR,
                 masked         = numberOfMaskedBases,
                 hits           = lengthRCoveredByHits,
                 hitsSA         = hitsSRSA$numberOfHits,
                 hitsSB         = hitsSRSB$numberOfHits,
                 hitsSAB        = numberOfDoubleHitPositions,
                 oldBkpBegin    = beginR + breakpointBegin,
                 oldBkpEnd      = beginR + breakpointEnd,
                 oldBkpLength   = oldBreakpointLength,
                 status         = ok,
                 left           = coefLeft,
                 right          = coefRight,
                 nonMaskedLeft  = coefLeftNonMasked,
                 nonMaskedRight = coefRightNonMasked,
                 bkpBegin       = beginR + newBreakpointBegin,
                 bkpEnd         = beginR + newBreakpointEnd,
                 bkpLength      = newBreakpointLength,
                 bkpNonMasked   = nonMaskedBreakpointLength,
                 bkpHits        = coveredByHitsBreakpointLength,
                 bkpDoubleHits  = numberOfDoubleHitBreakpointPositions,
                 RSS            = signif(segmentation$RSS, 3),
                 cbhBkpBegin    = coveredByHitsBreakpointBegin,
                 cbhBkpEnd      = coveredByHitsBreakpointEnd,
                 cbhLeft        = coefLeftCoveredByHits,
                 cbhRight       = coefRightCoveredByHits) 
    return (segmentation)
    
  } else {

    # If status == 1, it means that it is a valid breakpoint that have no
    # hits on the alignments SR vs SA and SR vs SB
    if (status == 1) {
      status = -1
    }
    
    # Calculate the old breakpoint length to add into the final dataframe
    oldBreakpointLength = breakpointEnd - breakpointBegin

    segmentation =
      data.frame(id             = breakpointId,
                 chrR           = chrR,
                 srLength       = lengthR,
                 masked         = numberOfMaskedBases,
                 hits           = 0,
                 hitsSA         = 0,
                 hitsSB         = 0,
                 hitsSAB        = 0,
                 oldBkpBegin    = beginR + breakpointBegin,
                 oldBkpEnd      = beginR + breakpointEnd,
                 oldBkpLength   = oldBreakpointLength,
                 status         = status,
                 left           = -1,
                 right          = -1,
                 nonMaskedLeft  = -1,
                 nonMaskedRight = -1,
                 bkpBegin       = beginR + breakpointBegin,
                 bkpEnd         = beginR + breakpointEnd,
                 bkpLength      = oldBreakpointLength,
                 bkpNonMasked   = 0,
                 bkpHits        = 0,
                 bkpDoubleHits  = 0,
                 RSS            = -1,
                 cbhBkpBegin    = -1,
                 cbhBkpEnd      = -1,
                 cbhLeft        = -1,
                 cbhRight       = -1) 
  } # if (status == 1 & (hitsSRSA$numberOfHits > 0 || hitsSRSB$numberOfHits > 0)) {...} else {...}
}

# #############################################################################
# FUNCTION transformIntervalsInPositions(intervalsTable)
#
# About the function:
#   This function receives a data frame that contains two columns: begin and
#   end positions of each alignment hit and creates a vector of positions that
#   are inside of these intervals. The positions are associated to the value 1
#   which indicates that in these positions, there is a hit on the alignment.
# #############################################################################
transformIntervalsInPositions = function(intervalsTable) {
  if (nrow(intervalsTable) > 0) {
    positionsVector = unlist(apply(intervalsTable, 1,
                                   function(x) seq(min(as.numeric(x[1]), as.numeric(x[2])),
                                                   max(as.numeric(x[1]), as.numeric(x[2])), by=1)),
                             use.names=FALSE)
    if (is.matrix(positionsVector)) {
      positionsVector = as.vector(positionsVector)
    }
  } else {
    positionsVector = vector()
  }
  return (unique(positionsVector))
}

# #############################################################################
# FUNCTION segmentMean(coveredByHitsScoreVector,
#                      coveredByHitsCumulatedIntervalLengths)
#
# About the function:
#  This function performs the segmentation. It receives a score vector of the
#  positions which are covered by hits and a vector that contains the cumulated
#  lengths of the intervals that are covered by hits.
#
#  The vector must be ordered and contain the 1 and n (length of y)
#
#  Here the model is more constrained, ie : on the first segment y = a with
#  a > 0, then y = 0 and on the third segment y = b with b < 0.
#
# #############################################################################
segmentMean = function(coveredByHitsScoreVector,
                        coveredByHitsCumulatedIntervalLengths) {

  # Get the number of positions of the vector coveredByHitsScoreVector
  numberOfPositions = length(coveredByHitsScoreVector)

  # Calculate the number of positions that have a value that is not zero
  sumall = sum(coveredByHitsScoreVector^2)

  # The cumulated interval lengths are the positions of the right extremities
  # of the intervals in the score vector (coveredByHitsScoreVector)
  endExtremities = coveredByHitsCumulatedIntervalLengths

  # The position 1 and all the position of the vector endExtremities
  # increased by one (except the last position) are the positions of the
  # left extremities of the intervals in the score vector
  # (coveredByHitsScoreVector)
  beginExtremities =
    c(1, endExtremities[1:length(endExtremities) - 1] + 1)
  
  # Calculate the cumulated sum of the score vector
  cumulatedScoreVectorSum = cumsum(coveredByHitsScoreVector)

  # Calculate the cumulated sum of the score vector (inverted)
  invCumulatedScoreVectorSum = cumsum(coveredByHitsScoreVector[numberOfPositions:1])[numberOfPositions:1]

  # Get the cumulated sum for the end extremity positions
  endExtremitiesCumulatedScore = cumulatedScoreVectorSum[endExtremities]

  # Get the cumulated sum for the begin extremity positions
  beginExtremitiesCumulatedScore = invCumulatedScoreVectorSum[beginExtremities]

  # Create a vector that hold the distance between a given begin extremity
  # and the end position
  distanceToTheEnd = numberOfPositions - beginExtremities + 1

  # Define the left and possibleLeft vectors
  left = 0
  possibleLeft = endExtremities[endExtremitiesCumulatedScore > 0]
  if (length(possibleLeft) > 0) {
    left = c(0, (endExtremitiesCumulatedScore^2/endExtremities)[endExtremitiesCumulatedScore > 0])
    possibleLeft = c(0, possibleLeft)
  } else {
    possibleLeft = 0
  }

  # Define the right and possibleRight vectors
  right = 0
  possibleRight = beginExtremities[beginExtremitiesCumulatedScore < 0]
  if (length(possibleRight) > 0) {
    right = c((beginExtremitiesCumulatedScore^2/distanceToTheEnd)[beginExtremitiesCumulatedScore < 0], 0)
    possibleRight = c(possibleRight, numberOfPositions + 1)
  } else{
    possibleRight = numberOfPositions + 1
  }

  # try if maximising independently left and right works :
  breakpointBegin = possibleLeft[tail(which.max(left), 1)]
  breakpointEnd = possibleRight[head(which.max(right), 1)]

  if (breakpointBegin < breakpointEnd) {
    # Just calculate the RSS value
    RSS = sumall - max(left) - max(right)
  } else{
    # Some aditional processing is necessary to correct the
    # breakpoint extremities
    nX = length(possibleLeft)

    tmpRSS = vector("numeric", nX)
    tmpPossibleRight = vector("numeric", nX)
    for (i in 1:nX) {
      maxRight = max(right[possibleRight > possibleLeft[i]])
      tmpRSS[i] = left[i] + maxRight
      tmpPossibleRight[i] =
        possibleRight[possibleRight > possibleLeft[i]][head(which.max(right[possibleRight > possibleLeft[i]]), 1)]
    }

    RSS = sumall - max(tmpRSS)
    ind = which.max(tmpRSS)
    breakpointBegin = possibleLeft[ind]
    breakpointEnd = tmpPossibleRight[ind]
    if (length(ind) > 1) {
      # If we have more than 1, we get the smallest one
      i = which.min(breakpointEnd - breakpointBegin)[1]
      breakpointBegin = breakpointBegin[i]
      breakpointEnd = breakpointEnd[i]
    }
  }
  return (list(RSS = RSS, begin = breakpointBegin, end = breakpointEnd))
}

# #############################################################################
# FUNCTION getTruePosition(positionOnSPrime, nonCoveredByHitsIntervals)
#
# About the function:
#  This function performs the correction of a position on S' to a position on
#  S (S' is the sequence S after the removal of the intervals that have no
#  alignment hits)
#
# Parameters:
#   positionOnSPrime = position on the sequence S'
#   nonCoveredByHitsIntervals = list of begin and end positions of the intervals
#                               (on the sequence S) that have no alignment hits
#
# #############################################################################
getTruePosition = function(positionOnSPrime, nonCoveredByHitsIntervals) {
  positionOnS = positionOnSPrime
  if (nrow(nonCoveredByHitsIntervals) > 0) {
    cumLength = c(0, cumsum(nonCoveredByHitsIntervals$end - nonCoveredByHitsIntervals$begin + 1))
    auxiliaryPosition = nonCoveredByHitsIntervals$begin - cumLength[-length(cumLength)]
    positionOnS = positionOnSPrime + cumLength[findInterval(positionOnSPrime, c(0, auxiliaryPosition))]
  }
  return (positionOnS)
}

# #############################################################################
# FUNCTION plotSegmentation(scoreTable, hitIntervalsSRSA, hitIntervalsSRSB,
#                           newBreakpointBegin, newBreakpointEnd, plotTitle,
#                           breakpointBegin, breakpointEnd,
#                           maskedIntervalsOnSR = NULL)
#
# About this function:
#  This function builds a plot of the segmentation that was performed on the
#  breakpoint region.
#
# #############################################################################
plotSegmentation = function(scoreTable, hitIntervalsSRSA, hitIntervalsSRSB,
                            newBreakpointBegin, newBreakpointEnd, plotTitle = "",
                            breakpointBegin = NULL, breakpointEnd = NULL,
                            maskedIntervalsOnSR = NULL) {

  hitSA = hitIntervalsSRSA
  hitSB = hitIntervalsSRSB

  y = cumsum(scoreTable$y)

  plot(c(min(scoreTable$x), max(scoreTable$x)),
       c(min(y), max(y)),
       type = "n",
       main = plotTitle,
       xlab = "Position on SR",
       ylab = "Cumulated score")

  minimumScore = min(y)

  # If we have the information about the masked sequence,
  # we draw light yellow rectangles to represent them
  if (!is.null(maskedIntervalsOnSR)) {
    if (nrow(maskedIntervalsOnSR) > 0) {
      maximumScore = max(y)
      rect(maskedIntervalsOnSR$begin, minimumScore,
           maskedIntervalsOnSR$end, maximumScore,
           col = "lightyellow2", border = NA)
    }
  }

  # Draw the score curve (black curve)
  lines(scoreTable$x, y)

  # Draw the original breakpoint extremities (vertical dark blue lines)
  if (!is.null(breakpointBegin)) {
    abline(v = breakpointBegin, col = 4)
  }
  if (!is.null(breakpointEnd)) {
    abline(v = breakpointEnd, col = 4)
  }

  # Draw the new breakpoint extremities (vertical red and green lines)
  abline(v = newBreakpointBegin, col = 2, lwd = 2)
  abline(v = newBreakpointEnd, col = 3, lwd = 2)

  # Draw the hit positions
  if (nrow(hitSA) > 0 & nrow(hitSB) > 0) {

    # Get positions that appears on the two alignments and draw light gray lines
    doubleHit = coveredBy(hitSB, hitSA)
    if (nrow(doubleHit) > 0) {
      segments(doubleHit$begin, minimumScore, doubleHit$end, minimumScore, col = 8, lwd = 10, lend = 1)
    }

    # Get positions that appears only on the alignmnets SR vs SA
    # and draw red lines
    uniqueA = notCoveredBy(hitSA, hitSB)
    if (nrow(uniqueA) > 0) {
      segments(uniqueA$begin, minimumScore, uniqueA$end, minimumScore, col = 2, lwd = 10, lend = 1)
    }

    # Get positions that appears only on the alignmnets SR vs SA
    # and draw green lines
    uniqueB = notCoveredBy(hitSB, hitSA)
    if (nrow(uniqueB) > 0) {
      segments(uniqueB$begin, minimumScore, uniqueB$end, minimumScore, col = 3, lwd = 10, lend = 1)
    }


  } else{
    # Draw red lines for the hits of the alignment SR vs SA
    if (nrow(hitSA) > 0) {
      segments(hitSA[,1], minimumScore, hitSA[,2], minimumScore, col = 2, lwd = 10, lend = 1)
    }
    # Draw green lines for the hits of the alignment SR vs SB
    if (nrow(hitSB) > 0) {
      segments(hitSB[,1], minimumScore, hitSB[,2], minimumScore, col = 3, lwd = 10, lend = 1)
    }
  }
}

# #############################################################################
# FUNCTION getHitsInformation(alignment, lengthSeqR)
#
# About this function:
#  This function receive a data frame with the information about the hits of
#  the alignment of the sequence SR against SA (or SB) and process it to
#  determine the intervals and positions that show a hit.
#
# #############################################################################
getHitsInformation = function(alignment, lengthSeqR) {

  # Data frame which holds the intervals that correspond to the hits of the
  # alignment between SR and the other sequence
  intervals = NULL

  # Vector which holds the positions that are inside of hit intervals of the
  # alignment between SR and the other sequence
  positions = NULL

  # Number of positions on the sequence SR that show a hit with the other
  # sequence
  numberOfHitPositions = 0

  if (nrow(alignment) == 0) {
    # Create an empty data frame of hits
    intervals = data.frame()

    # Create an vector that has lenght lengthSeqR full of zeros (it
    # indicates that there are no position on the sequence SR which
    # have hit with the other sequence)
    positions = rep(0, lengthSeqR)
  } else {
    # Create a vector of positions that are inside of the hits
    intervals = unionIntervals(alignment[,c(2,4)], T)

    # Transform the intervals and positions
    hitPositions = transformIntervalsInPositions(intervals)

    # Update the number of hit positions
    numberOfHitPositions = length(hitPositions)

    # Create a vector that has lenght "lengthSeqR" and have the value 0 for the
    # positions which show no hit and 1 for the positions that show a hit in
    # the alignment of the sequences SR and the other sequence
    positions = hist(hitPositions, breaks = seq(1, lengthSeqR + 1, by = 1),
                     plot = FALSE, right = F)$counts
  }

  return (list(intervals = intervals,
               positions = positions,
               numberOfHits = numberOfHitPositions))
}

# #############################################################################
# FUNCTION getCoveredByHitsCumulatedIntervalLengths(scoreVectorRuns)
#
# About this function:
#  This function get the score vector runs and return a vector of cumulated 
#  lengths of intervals that are covered by hits.
#
# #############################################################################
getCoveredByHitsCumulatedIntervalLengths = function(scoreVectorRuns) {

  # Get the vector of values
  v = scoreVectorRuns$values

  # Get the vector of lengths
  l = scoreVectorRuns$lengths

  # Add 0 to the last position of the vector to help computing all
  # cumulated intervals 
  v[length(v) + 1] = 0
  l[length(l) + 1] = 0

  # Create a data.frame with the two vectors
  data = data.frame(value=v, runlength=l)

  # Remove all lines that have value 2 (no hits with sequences A and B)
  # It means that just the intervals that are covered by hits will be
  # considered
  data = data[data$value != 2, ]

  # Compute the cumulated intervals lengths
  return (cumsum(data$runlength)[-length(data$runlength)])
}

# #############################################################################
# FUNCTION getNonCoveredByHitsIntervals(scoreVectorRuns)
#
# About this function:
#  This function gets the score vector and produces a data frame that have
#  all begin and end positions of the intervals that are non covered by hits.
#
# #############################################################################
getNonCoveredByHitsIntervals = function(scoreVectorRuns) {
  
  # Data frame that holds the begin and end positions of all intervals that
  # are not covered by any hit (consecutive elements that have the value 2 on
  # vector scoreVector)
  nonCoveredByHitsIntervals =
    data.frame(begin = cumsum(c(1, scoreVectorRuns$lengths))[c(scoreVectorRuns$values==2, FALSE)],
               end = cumsum(scoreVectorRuns$lengths)[scoreVectorRuns$values==2])

  return (nonCoveredByHitsIntervals)
}

# #############################################################################
# FUNCTION segmentMean(scoreVectorRuns, RSS)
#
# About the function:
#  This function performs the segmentation. It receives a score vector of the
#  positions which are covered by hits and a vector that contains the cumulated
#  lengths of the intervals that are covered by hits.
#
#  The vector must be ordered and contain the 1 and n (length of y)
#
#  Here the model is more constrained, ie : on the first segment y = a with
#  a > 0, then y = 0 and on the third segment y = b with b < 0.
#
# #############################################################################
statisticalTest = function(scoreVectorRuns, RSS) {


  # Get the vector of values
  v = scoreVectorRuns$values

  # Get the vector of lengths
  l = scoreVectorRuns$lengths

  # Create a data.frame with the two vectors
  data = data.frame(value=v, runlength=l)

  # Remove all lines that have value 2 (no hits with sequences A and B)
  # It means that just the intervals that are covered by hits will be
  # considered
  data = data[data$value != 2, ]

  n = 0
  
  for (i in seq(1:100)) {

    # Randomize the data
    randomizedData = randomizeScoreVector(data)

    # Calculate the segmentation over the randomized data
    segmentation = segmentMean(randomizedData$scoreVector, randomizedData$cumulatedIntervalLengths)

    # If the random data has a better RSS, we increase the counter.
    # If the counter reaches 5, it means that the original data cannot
    # be segmented
    if (RSS > segmentation$RSS) {
      n = n + 1
      if (n == 5) {
        return (0)
      }
    }
  }
  return (1)
}

# #############################################################################
# FUNCTION randomizeScoreVector(data)
#
# About the function:
#  This function receives a data frame that have the size and the value of each
#  interval and returns a list that has a randomized score vector and its
#  equivalent vector of cumulated interval lengths.
#
# #############################################################################
randomizeScoreVector = function(data) {

  # Get the number of intervals
  nIntervals = nrow(data)
  # Generate a random order of these intervals
  randomOrder = sample(seq(1:nIntervals))

  # Create the vector that will handle the new score vector and its equivalent
  # cumulated interval lengths
  newScoreVector = vector()
  newCumulatedIntervalLengthsVector = vector()

  # Auxiliary variables
  index = 1
  sum = 0

  while (index <= nIntervals) {

    # Get the selected interval (in the random order)
    selectedInterval = randomOrder[index]

    # Get the size and the value of this interval
    size = data$runlength[selectedInterval]
    value = data$value[selectedInterval]

    # Accumulate the size of the interval
    sum = sum + size

    # Add a interval of length "size" and value "value" in the new score vector
    newScoreVector = c(newScoreVector, seq(value, value, length.out = size))
    # Add the cumulated sum in the vector of cumulated interval lengths
    newCumulatedIntervalLengthsVector = c(newCumulatedIntervalLengthsVector, sum)

    index = index + 1
  }

  return (list(scoreVector = newScoreVector, cumulatedIntervalLengths = newCumulatedIntervalLengthsVector))
}
