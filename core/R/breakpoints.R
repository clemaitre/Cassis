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
# This file contains the definition of the function identifyBreakpoints
# Author: Claire Lemaitre
# Reviewed by: Christian Baudet
# #############################################################################

# #############################################################################
# FUNCTION identifyBreakpoints1(orthologyTable,
#                               limitR, limitAB,
#                               chromosomeBoundariesGO,
#                               extended, extendedAB,
#                               extendBeforeVerifyLength,
#                               minimumSequenceSize)
#
# About this function:
#
# 1 - This function receives an orthology table and identify all breakpoints
#
# 2 - This function receives the following parameters:
#     + orthologyTable - Table that contains the pairs of orthologs genes
#       between two genomes GR and GO
#     + limitR - Positive integer value that defines the maximum size for the
#       sequence R (Default value = 1000000000)
#     + limitAB - Positive integer value that defines the maximum size for the
#       sequences A and B (Default value = 1500000)
#     + chromosomeBoundariesGO - Table that have the information about the
#       begin and end position of each chromosome on the genome GO
#     + extended - Boolean parameter. If true, the genes A and B will be
#       included on the sequences R, A and B. If false, the sequences start/end
#       before/after the genes extremities (Default value = F)
#     + extendedAB - Boolean parameter. If true, the sequences A and B will
#       include the next/previous gene of the sequence. If false, the sequences
#       A and B start/finish after the next/previous gene.
#     + extendBeforeVerifyLength - Boolean parameter. If true, the method 
#       extend the sequences before verifying if they must be discarded,
#       otherwise, it verifies before applying the extension.
#     + minimumSequenceSize - Non negative integer parameter. Indicates the
#       minimum size that the sequences SR, SA and SB must have.
#
# 3 - The orthology table must have the following format (name and order):
#     g1      - Name of the gene in the genome GR
#     c1      - Chromosome where the gene g1 is located
#     inf1    - Start position of the gene g1 on the chromosome c1
#     sup1    - End position of the gene g1 on the chromosome c1
#     strand1 - Strand of the gene g1
#     g2      - Name of the gene in the genome GO
#     c2      - Chromosome where the gene g2 is located
#     inf2    - Start position of the gene g2 on the chromosome c2
#     sup2    - End position of the gene g2 on the chromosome c2
#     strand2 - Strand of the gene g2
#
# 4 - The chromosome boundaries table must have the following format
#     (name and order):
#     chr - Name of the chromosome
#     inf - Leftmost extremity of the chromosome
#     sup - Rightmost extremity of the chromosome
#
# 5 - The function will produce a data frame that has the following
#     columns:
#     id        - Breakpoint ID
#     type      - Type of the breakpoint: inter or intra
#     sRgeneA   - Name of the gene A on the sequence SR (genome GR)
#     sRgeneB   - Name of the gene B on the sequence SR (genome GR)
#     sRchr     - Chromosome of the genes A and B (genome GR)
#     sRstrandA - Strand of the gene A (genome GR)
#     sRstrandB - Strand of the gene B (genome GR)
#     sRinf     - Inferior boundary of the sequence SR
#     sRsup     - Superior boundary of the sequence SR
#     sAgene    - Name of the gene A on the sequence SA (genome GO)
#     sBgene    - Name of the gene B on the sequence SB (genome GO)
#     sAchr     - Chromosome of the gene A (genome GO)
#     sBchr     - Chromosome of the gene B (genome GO)
#     sAstrand  - Strand of the gene A (genome GO)
#     sBstrand  - Strand of the gene B (genome GO)
#     sAinf     - Inferior boundary of the sequence SA
#     sAsup     - Superior boundary of the sequence SA
#     sBinf     - Inferior boundary of the sequence SB
#     sBsup     - Superior boundary of the sequence SB
#     bkpBegin  - Relative position of the breakpoint begin (related to sRinf)
#     bkpEnd    - Relative position of the breakpoint end (related to sRinf)
#     status    - Status of the breakpoint
#
# 6 - Breakpoint status values:
#    1 : Valid breakpoint
#   -2 : Sequence SR smaller than the allowed limit
#   -3 : Sequence SA smaller than the allowed limit
#   -4 : Sequence SB smaller than the allowed limit
#   -5 : Sequence SA and SB smaller than the allowed limit
#   -6 : Sequence SR bigger than the allowed limit
#
# #############################################################################
identifyBreakpoints1 = function(orthologyTable,
                                limitR=1000000000, limitAB=3000000,
                                chromosomeBoundariesGO=NULL,
                                extended=F, extendedAB=F,
                                extendBeforeVerifyLength=T,
                                minimumSequenceSize=1) {

  # Avoid negative value for the parameter minimumSequenceSize
  if (minimumSequenceSize <= 0) {
    minimumSequenceSize = 1
  }

  # Sort the orthology table by the chromosome and the begining position
  # of the gene on the genome GR
  sortedByC1Inf1 = orthologyTable[order(orthologyTable$c1, orthologyTable$inf1),]

  # Sort the orthology table by the chromosome and the begining position
  # of the gene on the genome GO
  sortedByC2Inf2 = orthologyTable[order(orthologyTable$c2, orthologyTable$inf2),]

  # Create a vector that says, for each line of the orthology table
  # sortedByC1Inf1, the rank of the row when we sort the table by the
  # chromosome and the begining of the gene on the genome GO
  per = order(order(sortedByC1Inf1$c2, sortedByC1Inf1$inf2))

  # Create a vector that says, for each line of the orthology table
  # sortedByC1Inf1, the signed rank of the row when we sort the table by
  # the chromosome and the begining of the gene on the genome GO. The
  # rank value will be positive if the genes on both genomes have the same
  # direction. Otherwise it will be negative
  perS = ifelse(sortedByC1Inf1$strand1==sortedByC1Inf1$strand2, per, -per)

  # Create a data frame that combines the row i with the row i+1 of the data
  # frames sortedByC1Inf1 and perS (row 1 with row 2, row 2 with row 3, etc)
  # The dataframe will have the columns:
  # g1, c1, inf1, sup1, strand1, g2, c2, inf2, sup2, strand2 (row i of sortedByC1Inf1)
  # g1.1, c1.1, inf1.1, sup1.1, strand1.1, g2.1, c2.1, inf2.1, sup2.1, strand2.1 (row i+1 of sortedByC1Inf1)
  # perS..length.perS.. (row i from perS)
  # perS..1. (row i+1 from perS)
  auxDataframe = data.frame(sortedByC1Inf1[-nrow(sortedByC1Inf1),], sortedByC1Inf1[-1,], perS[-length(perS)], perS[-1])

  # Apply the defined function on each row of the dataframe to determine
  # the boundaries of the sequences A and B
  sequenceAB = apply(auxDataframe, 1, function(x) {

    # Determine if the breakpoint status
    #    1 : Valid breakpoint
    #   -2 : Sequence SR smaller than the allowed limit
    #   -3 : Sequence SA smaller than the allowed limit
    #   -4 : Sequence SB smaller than the allowed limit
    #   -5 : Sequence SA and SB smaller than the allowed limit
    #   -6 : Sequence SR bigger than the allowed limit
    #  -10 : Consecutive blocks of the genome GR are in different chromossomes
    #  -11 : We do not have a rearrangement
    ok = 1

    # Inferior boundary of the sequence R
    sequenceRinf = 0
    # Superior boundary of the sequence R
    sequenceRsup = 0
    # Inferior boundary of the sequence A
    sequenceAinf = 0
    # Superior boundary of the sequence A
    sequenceAsup = 0
    # Inferior boundary of the sequence B
    sequenceBinf = 0
    # Superior boundary of the sequence B
    sequenceBsup = 0

    # Non extended inferior boundary of the sequence R
    nonExtendedSequenceRinf = 0
    # Non extended superior boundary of the sequence R
    nonExtendedSequenceRsup = 0
    # Non extended inferior boundary of the sequence A
    nonExtendedSequenceAinf = 0
    # Non extended superior boundary of the sequence A
    nonExtendedSequenceAsup = 0
    # Non extended inferior boundary of the sequence B
    nonExtendedSequenceBinf = 0
    # Non extended superior boundary of the sequence B
    nonExtendedSequenceBsup = 0

    if ( x[2] != x[12] ) {
      # Two consecutives blocks on the genome GR are on different
      # chromosomes (we do not have a rearrangement here)
      ok = -10
    } else {
      # Two consecutives blocks on the genome GR are on the same
      # chromosome

      # Strand of the gene g1 on the first pair (g1,g2)
      pair1StrandG1 = as.numeric(x[5])
      # Strand of the gene g1 on the second pair (g1,g2)
      pair2StrandG1 = as.numeric(x[15])
      # Signed position of the gene g2 on the first pair (g1,g2)
      pair1SignedPositionG2 = as.numeric(x[21])
      # Signed position of the gene g2 on the second pair (g1,g2)
      pair2SignedPositionG2 = as.numeric(x[22])
      # Strand of the gene g2 on the first pair (g1,g2)
      pair1StrandG2 = as.numeric(x[10])
      # Strand of the gene g2 on the second pair (g1,g2)
      pair2StrandG2 = as.numeric(x[20])
      # Chromosome where is located the gene g2 of the first pair (g1,g2)
      pair1ChromosomeG2 = x[7]
      # Chromosome where is located the gene g2 of the second pair (g1,g2)
      pair2ChromosomeG2 = x[17]

      # Start position of the block g1 on the first pair
      pair1InfG1 = as.numeric(x[3])
      # End position of the block g1 on the first pair
      pair1SupG1 = as.numeric(x[4])
      # Start position of the block g2 on the second pair
      pair2InfG1 = as.numeric(x[13])
      # End position of the block g2 on the second pair
      pair2SupG1 = as.numeric(x[14])

      # Start position of the block g1 on the first pair
      pair1InfG2 = as.numeric(x[8])
      # End position of the block g1 on the first pair
      pair1SupG2 = as.numeric(x[9])
      # Start position of the block g2 on the second pair
      pair2InfG2 = as.numeric(x[18])
      # End position of the block g2 on the second pair
      pair2SupG2 = as.numeric(x[19])

      if (extended) {
        # If extended, the sequence SR includes the genes g1 and g2
        sequenceRinf = pair1InfG1
        sequenceRsup = pair2SupG1
      } else {
        # If not extended, the sequence SR includes just the region
        # between g1 and g2
        sequenceRinf = pair1SupG1
        sequenceRsup = pair2InfG1
      }
      nonExtendedSequenceRinf = pair1SupG1
      nonExtendedSequenceRsup = pair2InfG1

      # Length of the sequence R
      sequenceRlength = sequenceRsup - sequenceRinf
      nonExtendedSequenceRlength = nonExtendedSequenceRsup - nonExtendedSequenceRinf

      # Verify if there is overlapping on the sequence SR
      if ((extendBeforeVerifyLength == T & sequenceRlength < minimumSequenceSize) |
          (extendBeforeVerifyLength == F & nonExtendedSequenceRlength < minimumSequenceSize)) {
        ok = -2
      } else {
        
        if (pair1ChromosomeG2 != pair2ChromosomeG2 |
            pair2SignedPositionG2 != pair1SignedPositionG2 + 1) {
          # We have a rearrangement in this case

          if (sequenceRlength <= limitR) {
            # The region size agrees with the imposed limit

            # Sum of the lengths of the sequences A e B
            sequencesABlength = ifelse(sequenceRlength > limitAB, sequenceRlength, limitAB)

            # Get the unsigned positions of the gene g2 on the
            # first and second pairs
            pair1UnsignedPositionG2 = abs(pair1SignedPositionG2)
            pair2UnsignedPositionG2 = abs(pair2SignedPositionG2)

            # Define the sequence A
            if (pair1StrandG1 == pair1StrandG2) {
              # In this cases, the genes g1 and g2 are on the same strand (on
              # the first pair)

              # Inferior extremity of the sequence A
              if (extended) {
                # Include the gene on the sequence A
                sequenceAinf = pair1InfG2
              } else {
                # Include just the sequence after the gene
                sequenceAinf = pair1SupG2
              }
              # Non extended inferior extremity of the sequence A
              nonExtendedSequenceAinf = pair1SupG2
              
              # Calculate the end of the sequence A
              sequenceAend = 0
              if ((pair1UnsignedPositionG2 == dim(sortedByC1Inf1)[1]) ||
                  (sortedByC2Inf2[pair1UnsignedPositionG2 + 1, 7] != pair1ChromosomeG2)) {
                # In this case, (pair1UnsignedPositionG2) is the last gene of
                # the chromosome (pair1ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence A
                  # ends at the position: sequenceAinf + sequencesABlength
                  sequenceAend = sequenceAinf + sequencesABlength
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # A ends at the end of the chromosome
                  sequenceAend = getChromosomeEnd(pair1ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended superior extremity of the sequence A
                nonExtendedSequenceAsup = sequenceAend

              } else {
                # In this case, (pair1UnsignedPositionG2) is not the last gene
                # of the chromosome pair1ChromosomeG2. The end of the
                # sequence A is on the begining of the next gene
                if (extendedAB) {
                  # In this case, we include the next gene
                  sequenceAend = sortedByC2Inf2[pair1UnsignedPositionG2 + 1, 9]
                } else {
                  # In this case, we stop the sequence A before the next gene
                  sequenceAend = sortedByC2Inf2[pair1UnsignedPositionG2 + 1, 8]
                }
                # Non extended superior extremity of the sequence A
                nonExtendedSequenceAsup = sortedByC2Inf2[pair1UnsignedPositionG2 + 1, 8]
                
              }

              # If the length of the sequence A is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the end position of the sequence A.
              if (sequenceAend - sequenceAinf > sequencesABlength) {
                sequenceAsup = sequenceAinf + sequencesABlength
              } else {
                sequenceAsup = sequenceAend
              }

            } else {
              # In this cases, the genes g1 and g2 are on different strands (on
              # the first pair)

              # Superior extremity of the sequence A
              if (extended) {
                # Include the gene on the sequence A
                sequenceAsup = pair1SupG2
              } else {
                # Include just the sequence before the gene
                sequenceAsup = pair1InfG2
              }
              # Non extended superior extremity of the sequence A
              nonExtendedSequenceAsup = pair1InfG2
              
              # Calculate the begin of the sequence A
              sequenceAbegin = 0
              if ((pair1UnsignedPositionG2 - 1 == 0) ||
                  (sortedByC2Inf2[pair1UnsignedPositionG2 - 1, 7] != pair1ChromosomeG2) ) {
                # In this case, (pair1UnsignedPositionG2) is the first gene of the
                # chromosome (pair1ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence A
                  # begin at the position: sequenceAsup - sequencesABlength
                  # If the value is negative, the sequence start at the position 1
                  sequenceAbegin = ifelse(sequenceAsup - sequencesABlength < 1, 1, sequenceAsup - sequencesABlength)
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # A begins at the begin of the chromosome
                  sequenceAbegin = getChromosomeBegin(pair1ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended inferior extremity of the sequence A
                nonExtendedSequenceAinf = sequenceAbegin
                                
              } else {
                # In this case, (pair1UnsignedPositionG2) is not the first gene
                # of the chromosome (pair1ChromosomeG2). The begin of the
                # sequence A is on the end of the previous gene
                if (extendedAB) {
                  # In this case, we include the previous gene
                  sequenceAbegin = sortedByC2Inf2[pair1UnsignedPositionG2 - 1 , 8]
                } else {
                  # In this case, we start the sequence A after the previous gene
                  sequenceAbegin = sortedByC2Inf2[pair1UnsignedPositionG2 - 1 , 9]
                }
                # Non extended inferior extremity of the sequence A
                nonExtendedSequenceAinf = sortedByC2Inf2[pair1UnsignedPositionG2 - 1 , 9]

              }
              
              # If the length of the sequence A is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the begin position of the sequence A.
              if (sequenceAsup - sequenceAbegin > sequencesABlength) {
                sequenceAinf = sequenceAsup - sequencesABlength
              } else{
                sequenceAinf = sequenceAbegin
              }

            }

            # Defines the sequence B
            if (pair2StrandG1 == pair2StrandG2) {
              # In this cases, the genes g1 and g2 are on the same strand (on
              # the second pair)

              # Superior extremity of the sequence B
              if (extended) {
                # Include the gene on the sequence B
                sequenceBsup = pair2SupG2
              } else {
                # Include just the sequence before the gene
                sequenceBsup = pair2InfG2
              }
              # Non extended superior extremity of the sequence B
              nonExtendedSequenceBsup = pair2InfG2

              # Calculate the begin of the sequence B
              sequenceBbegin = 0
              if ((pair2UnsignedPositionG2 - 1 == 0) ||
                  (sortedByC2Inf2[pair2UnsignedPositionG2 - 1, 7] != pair2ChromosomeG2) ) {
                # In this case, (pair2UnsignedPositionG2) is the first gene of
                # the chromosome (pair2ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence B
                  # begin at the position: sequenceBsup - sequencesABlength
                  # If the value is negative, the sequence start at the position 1
                  sequenceBbegin = ifelse(sequenceBsup - sequencesABlength < 1, 1, sequenceBsup - sequencesABlength)
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # B begins at the begin of the chromosome
                  sequenceBbegin = getChromosomeBegin(pair2ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended inferior extremity of the sequence B
                nonExtendedSequenceBinf = sequenceBbegin
                
              } else {
                # In this case, (pair2UnsignedPositionG2) is not the first
                # gene of the chromosome (pair2ChromosomeG2). The begin of the
                # sequence B is on the end of the previous gene
                
                if (extendedAB) {
                  # In this case, we include the previous gene
                  sequenceBbegin = sortedByC2Inf2[pair2UnsignedPositionG2 - 1 , 8]
                } else {
                  # In this case, we start the sequence B after the previous gene
                  sequenceBbegin = sortedByC2Inf2[pair2UnsignedPositionG2 - 1 , 9]
                }
                # Non extended inferior extremity of the sequence B
                nonExtendedSequenceBinf = sortedByC2Inf2[pair2UnsignedPositionG2 - 1 , 9]
                
              }

              # If the length of the sequence B is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the begin position of the sequence B.
              if (sequenceBsup - sequenceBbegin > sequencesABlength) {
                sequenceBinf = sequenceBsup - sequencesABlength
              } else {
                sequenceBinf = sequenceBbegin
              }

            } else {
              # In this cases, the genes g1 and g2 are on different strands (on
              # the second pair)

              # Inferior extremity of the sequence B
              if (extended) {
                # Include the gene on the sequence B
                sequenceBinf = pair2InfG2
              } else {
                # Include just the sequence after the gene
                sequenceBinf = pair2SupG2
              }
              # Non extended inferior extremity of the sequence B
              nonExtendedSequenceBinf = pair2SupG2
              
              # Calculate the end of the sequence B
              sequenceBend = 0
              if ((pair2UnsignedPositionG2 == dim(sortedByC1Inf1)[1]) ||
                  (sortedByC2Inf2[pair2UnsignedPositionG2 + 1, 7] != pair2ChromosomeG2) ) {
                # In this case, (pair2UnsignedPositionG2) is the last gene of the
                # chromosome (pair2ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence B
                  # ends at the position: sequenceBinf + sequencesABlength
                  sequenceBend = sequenceBinf + sequencesABlength
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # B ends at the end of the chromosome
                  sequenceBend = getChromosomeEnd(pair2ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended superior extremity of the sequence B
                nonExtendedSequenceBsup = sequenceBend
                
              } else {
                # In this case, (pair2UnsignedPositionG2) is not the last gene
                # of the chromosome pair2ChromosomeG2. The end of the
                # sequence B is on the begining of the next gene

                if (extendedAB) {
                  # In this case, we include the next gene
                  sequenceBend = sortedByC2Inf2[pair2UnsignedPositionG2 + 1, 9]
                } else {
                  # In this case, we stop the sequence A before the next gene
                  sequenceBend = sortedByC2Inf2[pair2UnsignedPositionG2 + 1, 8]
                }
                # Non extended superior extremity of the sequence B
                nonExtendedSequenceBsup = sortedByC2Inf2[pair2UnsignedPositionG2 + 1, 8]
                
              }

              # If the length of the sequence B is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the end position of the sequence B.
              if (sequenceBend - sequenceBinf > sequencesABlength) {
                sequenceBsup = sequenceBinf + sequencesABlength
              } else {
                sequenceBsup = sequenceBend
              }
              
            }

            # Verify if the sequences SA and SB are valid sequences
            if (extendBeforeVerifyLength == T) {
              if (sequenceAsup - sequenceAinf < minimumSequenceSize &
                  sequenceBsup - sequenceBinf >= minimumSequenceSize) {
                ok = -3
              } else {
                if (sequenceBsup - sequenceBinf < minimumSequenceSize &
                    sequenceAsup - sequenceAinf >= minimumSequenceSize) {
                  ok = -4
                } else {
                  if (sequenceAsup - sequenceAinf < minimumSequenceSize &
                      sequenceBsup - sequenceBinf < minimumSequenceSize) {
                    ok = -5
                  }
                }
              }
            } else {
              if (nonExtendedSequenceAsup - nonExtendedSequenceAinf < minimumSequenceSize &
                  nonExtendedSequenceBsup - nonExtendedSequenceBinf >= minimumSequenceSize) {
                ok = -3
              } else {
                if (nonExtendedSequenceBsup - nonExtendedSequenceBinf < minimumSequenceSize &
                    nonExtendedSequenceAsup - nonExtendedSequenceAinf >= minimumSequenceSize) {
                  ok = -4
                } else {
                  if (nonExtendedSequenceAsup - nonExtendedSequenceAinf < minimumSequenceSize &
                      nonExtendedSequenceBsup - nonExtendedSequenceBinf < minimumSequenceSize) {
                    ok = -5
                  }
                }
              }
            }
            
          } else {
            # The sequence length does not agree with the imposed limit
            ok = -6
          } # if (sequenceRlength <= limitR) {...} else {...}

        } else {
          # We do not have a rearrangement in this case
          ok = -11
        } # if (pair1ChromosomeG2 != pair2ChromosomeG2 |
          #     pair2SignedPositionG2 != pair1SignedPositionG2 + 1) {...} else {...}

      } # if ((extendBeforeVerifyLength == T & sequenceRlength < minimumSequenceSize) |
        #     (extendBeforeVerifyLength == F & nonExtendedSequenceRlength < minimumSequenceSize)) {...} else {...}

    } # if ( x[2] != x[12] ) {...} else {...}

    # Return the flag that says if the interval is ok and the limits of the
    # sequences A, B and R
    return (c(ok, sequenceAinf, sequenceAsup, sequenceBinf, sequenceBsup, sequenceRinf, sequenceRsup))
  })

  # Get the inferior and superior positions of the sequences A and B (just
  # the rows that are classified with ok > -10)
  status       = sequenceAB[1, sequenceAB[1,] > -10]
  sequenceAinf = sequenceAB[2, sequenceAB[1,] > -10]
  sequenceAsup = sequenceAB[3, sequenceAB[1,] > -10]
  sequenceBinf = sequenceAB[4, sequenceAB[1,] > -10]
  sequenceBsup = sequenceAB[5, sequenceAB[1,] > -10]
  sequenceRinf = sequenceAB[6, sequenceAB[1,] > -10]
  sequenceRsup = sequenceAB[7, sequenceAB[1,] > -10]

  # Create a data frame with all breakpoint data of the rows that are
  # classified with ok > -10
  allBreakpointData = auxDataframe[sequenceAB[1,] > -10, ]

  # Calculate the begin and end positions of the breakpoint
  # (position relative to the begining of the sequence R)
  breakpointBegin = allBreakpointData[,4] - sequenceRinf
  breakpointEnd = allBreakpointData[,13] - sequenceRinf

  # Organize the data
  breakpoint = data.frame(allBreakpointData[,c(1,11,2,5,15)], sequenceRinf, sequenceRsup,
                          allBreakpointData[,c(6,16,7,17,10,20)], sequenceAinf, sequenceAsup,
                          sequenceBinf, sequenceBsup, status)

  # Compute the number of breakpoints
  numberOfBreakpoints = nrow(breakpoint)

  # Put an identifier for each row and add the information about
  # the begin and end of the breakpoin (relative positions)
  toReturn = data.frame(1:numberOfBreakpoints,
                        breakpoint,
                        breakpointBegin,
                        breakpointEnd)

  # Rename the columns
  colnames(toReturn) = c("id",
                         "sRgeneA", "sRgeneB", "sRchr", "sRstrandA", "sRstrandB", "sRinf", "sRsup",
                         "sAgene", "sBgene", "sAchr", "sBchr", "sAstrand", "sBstrand",
                         "sAinf", "sAsup", "sBinf", "sBsup", "status",
                         "bkpBegin", "bkpEnd")

  # Classify the rows
  toReturn$type = ifelse(toReturn$sAchr == toReturn$sBchr, "intra", "inter")

  # Reorganize the data
  toReturn = toReturn[,c("id", "type",
                         "sRgeneA", "sRgeneB", "sRchr", "sRstrandA", "sRstrandB", "sRinf", "sRsup",
                         "sAgene", "sBgene", "sAchr", "sBchr", "sAstrand", "sBstrand",
                         "sAinf", "sAsup", "sBinf", "sBsup",
                         "bkpBegin", "bkpEnd", "status")]
  return (toReturn)
}

# #############################################################################
# FUNCTION identifyBreakpoints2(orthologyTable,
#                               limitR, limitAB,
#                               chromosomeBoundariesGR,
#                               chromosomeBoundariesGO,
#                               extended, extendedAB,
#                               extendBeforeVerifyLength,
#                               minimunSequenceSize)
#
# About this function:
#
# 1 - This function receives an orthology table and identify all breakpoints
#
# 2 - This function receives the following parameters:
#     + orthologyTable - Table that contains the pairs of orthologs blocks
#       between the genomes GR and GO
#     + limitR - Positive integer value that defines the maximum size for the
#       sequence R (Default value = 1000000000)
#     + limitAB - Positive integer value that defines the maximum size for the
#       sequences A and B (Default value = 3000000)
#     + chromosomeBoundariesGR - Table that have the information about the begin
#       and end position of each chromosome on the genome GR
#     + chromosomeBoundariesGO - Table that have the information about the begin
#       and end position of each chromosome on the genome GO
#     + extended - Non negative integer parameter. If different of zero, the 
#       definition of the boundaries of the sequence SR, SA and, SB includes
#       a fragment of size "extended" of the blocks (Default value = 50000).
#     + extendedAB - Non negative integer parameter. If different of zero, the 
#       definition of the boundaries of the sequence SA and SB includes a 
#       fragment of size "extendedAB" of the previous/next block of the
#       concerned chromosome sequence (Default value = 50000).
#     + extendBeforeVerifyLength - Boolean parameter. If true, the method 
#       extend the sequences before verifying if they must be discarded,
#       otherwise, it verifies before applying the extension.
#     + minimumSequenceSize - Non negative integer parameter. Indicates the
#       minimum size that the sequences SR, SA and SB must have.
#
# 3 - The orthology table must have the following format (name and order):
#     g1      - Name of the block in the genome GR
#     c1      - Chromosome where the block g1 is located
#     inf1    - Start position of the block g1 on the chromosome c1
#     sup1    - End position of the block g1 on the chromosome c1
#     strand1 - Strand of the block g1
#     g2      - Name of the block in the genome GO
#     c2      - Chromosome where the block g2 is located
#     inf2    - Start position of the block g2 on the chromosome c2
#     sup2    - End position of the block g2 on the chromosome c2
#     strand2 - Strand of the block g2
#
# 4 - The chromosome boundaries table must have the following format
#     (name and order):
#     chr - Name of the chromosome
#     inf - Leftmost extremity of the chromosome
#     sup - Rightmost extremity of the chromosome
#
# 5 - The function will produce a data frame that has the following
#     columns:
#     id        - Breakpoint ID
#     type      - Type of the breakpoint: inter or intra
#     sRgeneA   - Name of the block A on the sequence sR (genome GR)
#     sRgeneB   - Name of the block B on the sequence sR (genome GR)
#     sRchr     - Chromosome of the blocks A and B (genome GR)
#     sRstrandA - Strand of the block A (genome GR)
#     sRstrandB - Strand of the block B (genome GR)
#     sRinf     - Inferior boundary of the sequence sR
#     sRsup     - Superior boundary of the sequence sR
#     sAgene    - Name of the block A on the sequence SA (genome GO)
#     sBgene    - Name of the block B on the sequence SB (genome GO)
#     sAchr     - Chromosome of the block A (genome GO)
#     sBchr     - Chromosome of the block B (genome GO)
#     sAstrand  - Strand of the block A (genome GO)
#     sBstrand  - Strand of the block B (genome GO)
#     sAinf     - Inferior boundary of the sequence SA
#     sAsup     - Superior boundary of the sequence SA
#     sBinf     - Inferior boundary of the sequence SB
#     sBsup     - Superior boundary of the sequence SB
#     bkpBegin  - Relative position of the breakpoint begin (related to sRinf)
#     bkpEnd    - Relative position of the breakpoint end (related to sRinf)
#     status    - Status of the breakpoint
#
# 6 - Breakpoint status values:
#    1 : Valid breakpoint
#   -2 : Sequence SR smaller than the allowed limit
#   -3 : Sequence SA smaller than the allowed limit
#   -4 : Sequence SB smaller than the allowed limit
#   -5 : Sequence SA and SB smaller than the allowed limit
#   -6 : Sequence SR bigger than the allowed limit
#
# #############################################################################
identifyBreakpoints2 = function(orthologyTable,
                                limitR=1000000000, limitAB=3000000,
                                chromosomeBoundariesGR=NULL,
                                chromosomeBoundariesGO=NULL,
                                extended = 50000, extendedAB = 50000,
                                extendBeforeVerifyLength=T,
                                minimumSequenceSize = 50000) {

  # Avoid negative value for the parameter minimumSequenceSize
  if (minimumSequenceSize <= 0) {
    minimumSequenceSize = 1
  }
  
  # Sort the orthology table by the chromosome and the begining position
  # of the block on the genome GR
  sortedByC1Inf1 = orthologyTable[order(orthologyTable$c1, orthologyTable$inf1),]
  
  # Sort the orthology table by the chromosome and the begining position
  # of the block on the genome GO
  sortedByC2Inf2 = orthologyTable[order(orthologyTable$c2, orthologyTable$inf2),]
  
  # Create a vector that says, for each line of the orthology table
  # sortedByC1Inf1, the rank of the row when we sort the table by the
  # chromosome and the begining of the block on the genome GO
  per = order(order(sortedByC1Inf1$c2, sortedByC1Inf1$inf2))

  # Create a vector that says, for each line of the orthology table
  # sortedByC1Inf1, the signed rank of the row when we sort the table by
  # the chromosome and the begining of the block on the genome GO
  # The rank value will be positive if the blocks of the two genomes
  # have the same direction and negative otherwise
  perS = ifelse(sortedByC1Inf1$strand1==sortedByC1Inf1$strand2, per, -per)

  # Create a data frame that combines the row i with the row i+1 of the data
  # frames sortedByC1Inf1 and perS (row 1 with row 2, row 2 with row 3, etc)
  # The dataframe will have the columns:
  # g1, c1, inf1, sup1, strand1, g2, c2, inf2, sup2, strand2 (row i of sortedByC1Inf1)
  # g1.1, c1.1, inf1.1, sup1.1, strand1.1, g2.1, c2.1, inf2.1, sup2.1, strand2.1 (row i+1 of sortedByC1Inf1)
  # perS..length.perS.. (row i from perS)
  # perS..1. (row i+1 from perS)
  auxDataframe = data.frame(sortedByC1Inf1[-nrow(sortedByC1Inf1),], sortedByC1Inf1[-1,], perS[-length(perS)], perS[-1])

  
  # Apply the defined function on each row of the dataframe to determine
  # the boundaries of the sequences A and B
  sequenceAB = apply(auxDataframe, 1, function(x) {

    # Determine if the breakpoint status
    #    1 : Valid breakpoint
    #   -2 : Sequence SR smaller than the allowed limit 
    #   -3 : Sequence SA smaller than the allowed limit 
    #   -4 : Sequence SB smaller than the allowed limit 
    #   -5 : Sequence SA and SB smaller than the allowed limit 
    #   -6 : Sequence SR bigger than the allowed limit
    #  -10 : Consecutive blocks of the genome GR are in different chromossomes
    #  -11 : We do not have a rearrangement
    ok = 1
    
    # Inferior boundary of the sequence R
    sequenceRinf = 0
    # Superior boundary of the sequence R
    sequenceRsup = 0
    # Inferior boundary of the sequence A
    sequenceAinf = 0
    # Superior boundary of the sequence A
    sequenceAsup = 0
    # Inferior boundary of the sequence B
    sequenceBinf = 0
    # Superior boundary of the sequence B
    sequenceBsup = 0

    # Non extended inferior boundary of the sequence R
    nonExtendedSequenceRinf = 0
    # Non extended superior boundary of the sequence R
    nonExtendedSequenceRsup = 0
    # Non extended inferior boundary of the sequence A
    nonExtendedSequenceAinf = 0
    # Non extended superior boundary of the sequence A
    nonExtendedSequenceAsup = 0
    # Non extended inferior boundary of the sequence B
    nonExtendedSequenceBinf = 0
    # Non extended superior boundary of the sequence B
    nonExtendedSequenceBsup = 0

    if ( x[2] != x[12] ) {
      # Two consecutives blocks from the genome GR are on different
      # chromosomes
      ok = -10
    } else {
      # Two consecutives blocks from the genome GR are on the same
      # chromosome
      
      # Strand of the block g1 on the first pair (g1,g2)
      pair1StrandG1 = as.numeric(x[5])
      # Strand of the block g1 on the second pair (g1,g2)
      pair2StrandG1 = as.numeric(x[15])
      # Signed position of the block g2 on the first pair (g1,g2)
      pair1SignedPositionG2 = as.numeric(x[21])
      # Signed position of the block g2 on the second pair (g1,g2)
      pair2SignedPositionG2 = as.numeric(x[22])
      # Strand of the block g2 on the first pair (g1,g2)
      pair1StrandG2 = as.numeric(x[10])
      # Strand of the block g2 on the second pair (g1,g2)
      pair2StrandG2 = as.numeric(x[20])
      # Chromosome where is located the block g2 of the first pair (g1,g2)
      pair1ChromosomeG2 = x[7]
      # Chromosome where is located the block g2 of the second pair (g1,g2)
      pair2ChromosomeG2 = x[17]

      # Start position of the block g1 on the first pair
      pair1InfG1 = as.numeric(x[3])
      # End position of the block g1 on the first pair
      pair1SupG1 = as.numeric(x[4])
      # Start position of the block g2 on the second pair
      pair2InfG1 = as.numeric(x[13])
      # End position of the block g2 on the second pair
      pair2SupG1 = as.numeric(x[14])

      # Start position of the block g1 on the first pair
      pair1InfG2 = as.numeric(x[8])
      # End position of the block g1 on the first pair
      pair1SupG2 = as.numeric(x[9])
      # Start position of the block g2 on the second pair
      pair2InfG2 = as.numeric(x[18])
      # End position of the block g2 on the second pair
      pair2SupG2 = as.numeric(x[19])

      if (extended > 0) {
        # If extended, the sequence SR includes some additional
        # fragments of the blocks g1s of the first and
        # second pairs
        sequenceRinf = pair1SupG1 - extended
        sequenceRsup = pair2InfG1 + extended
        sequenceRinf = validateBeginPosition(sequenceRinf, x[2], chromosomeBoundariesGR)
        sequenceRsup = validateEndPosition(sequenceRsup, x[2], chromosomeBoundariesGR)
      } else {
        # If not extended, the sequence SR includes just the region
        # between the blocks g1s of the first and
        # second pairs
        sequenceRinf = pair1SupG1
        sequenceRsup = pair2InfG1
      }
      nonExtendedSequenceRinf = pair1SupG1
      nonExtendedSequenceRsup = pair2InfG1

      # Length of the sequence R
      sequenceRlength = sequenceRsup - sequenceRinf
      nonExtendedSequenceRlength = nonExtendedSequenceRsup - nonExtendedSequenceRinf
      
      # Verify if there is overlapping on the sequence SR
      if ((extendBeforeVerifyLength == T & sequenceRlength < minimumSequenceSize) |
          (extendBeforeVerifyLength == F & nonExtendedSequenceRlength < minimumSequenceSize)) {
        ok = -2
      } else {
        
        if (pair1ChromosomeG2 != pair2ChromosomeG2 |
            pair2SignedPositionG2 != pair1SignedPositionG2 + 1) {
          # We have a rearrangement in this case

          if (sequenceRlength <= limitR) {
            # The sequence length agrees with the imposed limit

            # Sum of the lengths of the sequences A e B
            sequencesABlength = ifelse(sequenceRlength > limitAB, sequenceRlength, limitAB)

            # Get the unsigned positions of the block g2 on the
            # first and second pairs
            pair1UnsignedPositionG2 = abs(pair1SignedPositionG2)
            pair2UnsignedPositionG2 = abs(pair2SignedPositionG2)

            # Define the sequence A
            if (pair1StrandG1 == pair1StrandG2) {
              # In this cases, the blocks g1 and g2 are on the same
              # strand (on the first pair)

              # Inferior extremity of the sequence A
              if (extended > 0) {
                # Include some additional fragment of the block g2
                # of the first pair
                sequenceAinf = pair1SupG2 - extended
                sequenceAinf = validateBeginPosition(sequenceAinf, pair1ChromosomeG2, chromosomeBoundariesGO)
              } else {
                # Include just the sequence after the block g2
                # of the first pair
                sequenceAinf = pair1SupG2
              }
              # Non extended inferior extremity of the sequence A
              nonExtendedSequenceAinf = pair1SupG2

              # Calculate the end of the sequence A
              sequenceAend = 0
              if ((pair1UnsignedPositionG2 == dim(sortedByC1Inf1)[1]) ||
                  (sortedByC2Inf2[pair1UnsignedPositionG2 + 1, 7] != pair1ChromosomeG2)) {
                # In this case, (pair1UnsignedPositionG2) is the last block
                # of the chromosome (pair1ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence A
                  # ends at the position: sequenceAinf + sequencesABlength
                  sequenceAend = sequenceAinf + sequencesABlength
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # A ends at the end of the chromosome
                  sequenceAend = getChromosomeEnd(pair1ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended superior extremity of the sequence A
                nonExtendedSequenceAsup = sequenceAend
              } else {
                # In this case, (pair1UnsignedPositionG2) is not the last 
                # block of the chromosome pair1ChromosomeG2. The end of
                # the sequence A is on the begining of the next block on
                # the genome GO
              
                nextPair1InfG2 = sortedByC2Inf2[pair1UnsignedPositionG2 + 1, 8]
                if (extendedAB > 0) {
                  # In this case, we include some fragment of the next block
                  # of the chromossome c2 on the genome GO
                  sequenceAend = nextPair1InfG2 + extendedAB
                  sequenceAend = validateEndPosition(sequenceAend, pair1ChromosomeG2, chromosomeBoundariesGO)
                } else {
                  # In this case, we stop the sequence A before the next block
                  # of the chromossome c2 on the genome GO
                  sequenceAend = nextPair1InfG2
                }
                # Non extended superior extremity of the sequence A
                nonExtendedSequenceAsup = nextPair1InfG2                
              }

              # If the length of the sequence A is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the end position of the sequence A.
              if (sequenceAend - sequenceAinf > sequencesABlength) {
                sequenceAsup = sequenceAinf + sequencesABlength
              } else {
                sequenceAsup = sequenceAend
              }

            } else {
              # In this cases, the blocks g1 and g2 are on different
              # strands (on the first pair)

              # Superior extremity of the sequence A
              if (extended > 0) {
                # Include some additional fragment of the block g2
                # of the first pair
                sequenceAsup = pair1InfG2 + extended
                sequenceAsup = validateEndPosition(sequenceAsup, pair1ChromosomeG2, chromosomeBoundariesGO)
              } else {
                # Include just the sequence before the block g2
                # of the first pair
                sequenceAsup = pair1InfG2
              }
              # Non extended superior extremity of the sequence A
              nonExtendedSequenceAsup = pair1InfG2

              # Calculate the begin of the sequence A
              sequenceAbegin = 0
              if ((pair1UnsignedPositionG2 - 1 == 0) ||
                  (sortedByC2Inf2[pair1UnsignedPositionG2 - 1, 7] != pair1ChromosomeG2) ) {
                # In this case, (pair1UnsignedPositionG2) is the first block
                # of the chromosome (pair1ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence A
                  # begin at the position: sequenceAsup - sequencesABlength
                  # If the value is negative, the sequence start at the position 1
                  sequenceAbegin = ifelse(sequenceAsup - sequencesABlength < 1, 1, sequenceAsup - sequencesABlength)
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # A begins at the begin of the chromosome
                  sequenceAbegin = getChromosomeBegin(pair1ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended inferior extremity of the sequence A
                nonExtendedSequenceAinf = sequenceAbegin
              
              } else {
                # In this case, (pair1UnsignedPositionG2) is not the first
                # block of the chromosome (pair1ChromosomeG2). The begin
                # of the sequence A is on the end of the previous block
              
                previousPair1SupG2 = sortedByC2Inf2[pair1UnsignedPositionG2 - 1, 9]
                if (extendedAB > 0) {
                  # In this case, we include some fragment of the previous block
                  # of the chromossome c2 on the genome GO
                  sequenceAbegin = previousPair1SupG2 - extendedAB
                  sequenceAbegin = validateBeginPosition(sequenceAbegin, pair1ChromosomeG2, chromosomeBoundariesGO)
                } else {
                  # In this case, we start the sequence A after the previous block
                  # of the chromossome c2 on the genome GO
                  sequenceAbegin = previousPair1SupG2
                }
                # Non extended inferior extremity of the sequence A
                nonExtendedSequenceAinf = previousPair1SupG2

              }

              # If the length of the sequence A is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the begin position of the sequence A.
              if (sequenceAsup - sequenceAbegin > sequencesABlength) {
                sequenceAinf = sequenceAsup - sequencesABlength
              } else{
                sequenceAinf = sequenceAbegin
              }

            }

            # Defines the sequence B
            if (pair2StrandG1 == pair2StrandG2) {
              # In this cases, the blocks g1 and g2 are on the same
              # strand (on the second pair)

              # Superior extremity of the sequence B
              if (extended > 0) {
                # Include some additional fragment of the block g2
                # of the second pair
                sequenceBsup = pair2InfG2 + extended
                sequenceBsup = validateEndPosition(sequenceBsup, pair2ChromosomeG2, chromosomeBoundariesGO)
              } else {
                # Include just the sequence before the block g2
                # of the second pair
                sequenceBsup = pair2InfG2
              }
              # Non extended superior extremity of the sequence B
              nonExtendedSequenceBsup = pair2InfG2

              # Calculate the begin of the sequence B
              sequenceBbegin = 0
              if ((pair2UnsignedPositionG2 - 1 == 0) ||
                  (sortedByC2Inf2[pair2UnsignedPositionG2 - 1, 7] != pair2ChromosomeG2) ) {
                # In this case, (pair2UnsignedPositionG2) is the first block
                # of the chromosome (pair2ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence B
                  # begin at the position: sequenceBsup - sequencesABlength
                  # If the value is negative, the sequence start at the position 1
                  sequenceBbegin = ifelse(sequenceBsup - sequencesABlength < 1, 1, sequenceBsup - sequencesABlength)
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # B begins at the begin of the chromosome
                  sequenceBbegin = getChromosomeBegin(pair2ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended inferior extremity of the sequence B
                nonExtendedSequenceBinf = sequenceBbegin
              
              } else {
                # In this case, (pair2UnsignedPositionG2) is not the first
                # block of the chromosome (pair2ChromosomeG2). The begin 
                # of the sequence B is on the end of the previous block

                previousPair1SupG2 = sortedByC2Inf2[pair2UnsignedPositionG2 - 1 , 9]
                if (extendedAB > 0) {
                  # In this case, we include some fragment of the previous block
                  # of the chromossome c2 on the genome GO
                  sequenceBbegin = previousPair1SupG2 - extendedAB
                  sequenceBbegin = validateBeginPosition(sequenceBbegin, pair2ChromosomeG2, chromosomeBoundariesGO)
                } else {
                  # In this case, we start the sequence B after the previous block
                  # of the chromossome c2 on the genome GO
                  sequenceBbegin = previousPair1SupG2
                }
                # Non extended inferior extremity of the sequence B
                nonExtendedSequenceBinf = previousPair1SupG2
              
              }

              # If the length of the sequence B is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the begin position of the sequence B.
              if (sequenceBsup - sequenceBbegin > sequencesABlength) {
                sequenceBinf = sequenceBsup - sequencesABlength
              } else {
                sequenceBinf = sequenceBbegin
              }

            } else {
              # In this cases, the blocks g1 and g2 are on different strands (on
              # the second pair)

              # Inferior extremity of the sequence B
              if (extended > 0) {
                # Include some additional fragment of the block g2
                # of the second pair
                sequenceBinf = pair2SupG2 - extended
                sequenceBinf = validateBeginPosition(sequenceBinf, pair2ChromosomeG2, chromosomeBoundariesGO)
              } else {
                # Include just the sequence after the block g2
                # of the second pair
                sequenceBinf = pair2SupG2
              }
              # Non extended inferior extremity of the sequence B
              nonExtendedSequenceBinf = pair2SupG2

              # Calculate the end of the sequence B
              sequenceBend = 0
              if ((pair2UnsignedPositionG2 == dim(sortedByC1Inf1)[1]) ||
                  (sortedByC2Inf2[pair2UnsignedPositionG2 + 1, 7] != pair2ChromosomeG2) ) {
                # In this case, (pair2UnsignedPositionG2) is the last block of the
                # chromosome (pair2ChromosomeG2)

                if (is.null(chromosomeBoundariesGO)) {
                  # If the table of chromosome length is NULL, the sequence B
                  # ends at the position: sequenceBinf + sequencesABlength
                  sequenceBend = sequenceBinf + sequencesABlength
                } else {
                  # If the table of chromosome length is not NULL, the sequence
                  # B ends at the end of the chromosome
                  sequenceBend = getChromosomeEnd(pair2ChromosomeG2, chromosomeBoundariesGO)
                }
                # Non extended inferior extremity of the sequence B
                nonExtendedSequenceBsup = sequenceBend
              
              } else {
                # In this case, (pair2UnsignedPositionG2) is not the last 
                # block of the chromosome pair2ChromosomeG2. The end
                # of the sequence B is on the begining of the next block

                nextPair1InfG2 = sortedByC2Inf2[pair2UnsignedPositionG2 + 1, 8]
                if (extendedAB > 0) {
                  # In this case, we include some fragment of the next block
                  # of the chromossome c2 on the genome GO
                  sequenceBend = nextPair1InfG2 + extendedAB
                  sequenceBend = validateEndPosition(sequenceBend, pair2ChromosomeG2, chromosomeBoundariesGO)
                } else {
                  # In this case, we stop the sequence A before the next block
                  # of the chromossome c2 on the genome GO
                  sequenceBend = nextPair1InfG2
                }
                # Non extended inferior extremity of the sequence B
                nonExtendedSequenceBsup = nextPair1InfG2

              }

              # If the length of the sequence B is bigger than the sum of the
              # expected sum of the lengths of the sequences A and B, correct
              # the end position of the sequence B.
              if (sequenceBend - sequenceBinf > sequencesABlength) {
                sequenceBsup = sequenceBinf + sequencesABlength
              } else {
                sequenceBsup = sequenceBend
              }
            }

            # Verify if the sequences SA and SB are valid sequences
            if (extendBeforeVerifyLength == T) {
              if (sequenceAsup - sequenceAinf < minimumSequenceSize &
                  sequenceBsup - sequenceBinf >= minimumSequenceSize) {
                ok = -3
              } else {
                if (sequenceBsup - sequenceBinf < minimumSequenceSize &
                    sequenceAsup - sequenceAinf >= minimumSequenceSize) {
                  ok = -4
                } else {
                  if (sequenceAsup - sequenceAinf < minimumSequenceSize &
                      sequenceBsup - sequenceBinf < minimumSequenceSize) {
                    ok = -5
                  }
                }
              }
            } else {
              if (nonExtendedSequenceAsup - nonExtendedSequenceAinf < minimumSequenceSize &
                  nonExtendedSequenceBsup - nonExtendedSequenceBinf >= minimumSequenceSize) {
                ok = -3
              } else {
                if (nonExtendedSequenceBsup - nonExtendedSequenceBinf < minimumSequenceSize &
                    nonExtendedSequenceAsup - nonExtendedSequenceAinf >= minimumSequenceSize) {
                  ok = -4
                } else {
                  if (nonExtendedSequenceAsup - nonExtendedSequenceAinf < minimumSequenceSize &
                      nonExtendedSequenceBsup - nonExtendedSequenceBinf < minimumSequenceSize) {
                    ok = -5
                  }
                }
              }
            }
            
          } else {
            # The sequence length does not agree with the imposed limit            
            ok = -6
          } # if (sequenceRlength <= limitR) {...} else {...}
          
        } else {
          # We do not have a rearrangement in this case
          ok = -11
        } #  if (pair1ChromosomeG2 != pair2ChromosomeG2 |
          #  pair2SignedPositionG2 != pair1SignedPositionG2 + 1) {...} else {...}
        
      } # if ((extendBeforeVerifyLength == T & sequenceRlength < minimumSequenceSize) |
        #     (extendBeforeVerifyLength == F & nonExtendedSequenceRlength < minimumSequenceSize)) {...} else {...}

    } # if ( x[2] != x[12] ) {...} else {...}

    # Return the flag that says if the interval is ok and the limits of the
    # sequences A and B
    return (c(ok, sequenceAinf, sequenceAsup, sequenceBinf, sequenceBsup, sequenceRinf, sequenceRsup))
    
  })

  # Get the inferior and superior positions of the sequences A and B (just
  # the rows that are classified with ok > -10)
  status       = sequenceAB[1, sequenceAB[1,] > -10]
  sequenceAinf = sequenceAB[2, sequenceAB[1,] > -10]
  sequenceAsup = sequenceAB[3, sequenceAB[1,] > -10]
  sequenceBinf = sequenceAB[4, sequenceAB[1,] > -10]
  sequenceBsup = sequenceAB[5, sequenceAB[1,] > -10]
  sequenceRinf = sequenceAB[6, sequenceAB[1,] > -10]
  sequenceRsup = sequenceAB[7, sequenceAB[1,] > -10]
  
  # Create a data frame with all breakpoint data of the rows that are
  # classified with ok > -10
  allBreakpointData = auxDataframe[sequenceAB[1,] > -10, ]

  # Calculate the begin and end positions of the breakpoint
  # (position relative to the begining of the sequence R)
  breakpointBegin = allBreakpointData[,4] - sequenceRinf
  breakpointEnd = allBreakpointData[,13] - sequenceRinf

  # Organize the data
  breakpoint = data.frame(allBreakpointData[,c(1,11,2,5,15)], sequenceRinf, sequenceRsup,
                          allBreakpointData[,c(6,16,7,17,10,20)], sequenceAinf, sequenceAsup,
                          sequenceBinf, sequenceBsup, status)

  # Compute the number of breakpoints
  numberOfBreakpoints = nrow(breakpoint)

  # Put an identifier for each row and add the information about
  # the begin and end of the breakpoin (relative positions)
  toReturn = data.frame(1:numberOfBreakpoints,
                        breakpoint,
                        breakpointBegin,
                        breakpointEnd)

  # Rename the columns
  colnames(toReturn) = c("id",
                         "sRgeneA", "sRgeneB", "sRchr", "sRstrandA", "sRstrandB", "sRinf", "sRsup",
                         "sAgene", "sBgene", "sAchr", "sBchr", "sAstrand", "sBstrand",
                         "sAinf", "sAsup", "sBinf", "sBsup", "status",
                         "bkpBegin", "bkpEnd")

  # Classify the rows
  toReturn$type = ifelse(toReturn$sAchr == toReturn$sBchr, "intra", "inter")

  # Reorganize the data
  toReturn = toReturn[,c("id", "type",
                         "sRgeneA", "sRgeneB", "sRchr", "sRstrandA", "sRstrandB", "sRinf", "sRsup",
                         "sAgene", "sBgene", "sAchr", "sBchr", "sAstrand", "sBstrand",
                         "sAinf", "sAsup", "sBinf", "sBsup",
                         "bkpBegin", "bkpEnd", "status")]
  return (toReturn)
}
