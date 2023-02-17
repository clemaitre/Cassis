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
# This file contains the definition of the function filterOverlappingGenes
#
# Author: Claire Lemaitre
# Reviewed by: Christian Baudet
# #############################################################################

# #############################################################################
# FUNCTION filterOverlappingGenes(orthologyTable, logFileName)
#
# About the function:
# 1 - This function search for overlapping genes.
#
# 2 - This function receives a table of ortholog genes which has the columns:
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
# 3 - The columns of the table must be named and ordered as described above.
#
# 4 - The input table can not have duplication of genes and must have
#     only one to one relationships between the genes of the two species.
#
# 5 - If the overlapping genes, which are in the same interval, are in
#     the same order and have the same direction in the two genomes,
#     they are merged to represent just 1 gene (this new gene will
#     receive name m+name of the first gene).
#
# 6 - If the overlapping genes do not match the above criteria, we
#     discard them.
#
# 7 - The function returns a list of tables:
#     + nonOverlappingGenesAndMergedIntervals = this dataframe contains the
#       list of non overlapping genes and the list of merged intervals. The
#       format of this dataframe is the same of the input dataframe.
#     + nonOverlappingGenes = this dataframe contains the list of non
#       overlapping genes. The format of this dataframe is the same of
#       the input dataframe.
#     + mergedIntervals = this dataframe contains the list of merged intervals.
#       The format of this dataframe is the same of the input dataframe.
#     + statistics = vector with statistics values about the intervals:
#       - numberOfOverlappingIntervalsG1 = number of overlapping intervals
#         that were found on the first species
#       - numberOfGenesInsideOfOverlappingIntervalsG1 = number of genes which
#         are inside of overlapping intervals that were found on the first
#         species
#       - numberOfMergedIntervalsG1 = number of merged intervals on the first
#         species
#       - numberOfGenesInsideOfMergedIntervalsG1 = number of genes which are
#         inside of merged intervals on the first species
#       - numberOfNonMergedIntervalsG1 = number of overlapping intervals on
#         the first species that could not be merged
#       - numberOfGenesInsideOfNonMergedIntervalsG1 = number of genes which
#         are inside of overlapping intervals on the first species that
#         could not be merged
#       - numberOfOverlappingIntervalsG2 = number of overlapping intervals
#         that were found on the second species
#       - numberOfGenesInsideOfOverlappingIntervalsG2 = number of genes which
#         are inside of overlapping intervals that were found on the second
#         species
#       - numberOfMergedIntervalsG2 = number of merged intervals on the second
#         species
#       - numberOfGenesInsideOfMergedIntervalsG2 = number of genes which are
#         inside of merged intervals on the second species
#       - numberOfNonMergedIntervalsG2 = number of overlapping intervals on
#         the second species that could not be merged
#       - numberOfGenesInsideOfNonMergedIntervalsG2 = number of genes which
#         are inside of overlapping intervals on the second species that
#         could not be merged
#     + genesList = a list of that have three list of genes:
#       - genesOnIntervals - List of genes that are inside of intervals
#       - genesOnMergedIntervals = List of genes that are inside of merged
#         intervals
#       - discardedGenes = List of genes that are inside of intervals which
#         could not be merged.
#
# 8 - If the parameter logFileName is not NULL, the function will register
#     some information about the steps it is performing in a file that has
#     the name which is specified by the parameter
#
# #############################################################################
filterOverlappingGenes = function(orthologyTable, logFileName=NULL) {

  # Warning! The calculus over the second species are performed after we
  #          had discarded the merged genes on the first step

  # Create the auxiliary vectors that will handle the merged intervals
  merge.g1 = vector()
  merge.c1 = vector()
  merge.inf1 = vector()
  merge.sup1 = vector()
  merge.strand1 = vector()
  merge.g2 = vector()
  merge.c2 = vector()
  merge.inf2 = vector()
  merge.sup2 = vector()
  merge.strand2 = vector()

  # This vector keeps the list of genes that are inside of intervals
  genesOnIntervals = vector()
  # This vector keeps the list of genes that are merged
  genesOnMergedIntervals = vector()

  # Number of overlapping intervals on the first genome
  numberOfOverlappingIntervalsG1 = 0
  # Number of gene inside of overlapping intervals on the first genome
  numberOfGenesInsideOfOverlappingIntervalsG1 = 0
  # Number of merged intervals on the first genome
  numberOfMergedIntervalsG1 = 0
  # Number of genes inside of merged intervals on the first genome
  numberOfGenesInsideOfMergedIntervalsG1 = 0
  # Number of intervals that could not be merged on the first genome
  numberOfNonMergedIntervalsG1 = 0
  # Number of genes inside of intervals that could not be merged on the first genome
  numberOfGenesInsideOfNonMergedIntervalsG1 = 0

  # Number of overlapping intervals on the second genome
  numberOfOverlappingIntervalsG2 = 0
  # Number of gene inside of overlapping intervals on the second genome
  numberOfGenesInsideOfOverlappingIntervalsG2 = 0
  # Number of merged intervals on the second genome
  numberOfMergedIntervalsG2 = 0
  # Number of genes inside of merged intervals on the second genome
  numberOfGenesInsideOfMergedIntervalsG2 = 0
  # Number of intervals that could not be merged on the second genome
  numberOfNonMergedIntervalsG2 = 0
  # Number of genes inside of intervals that could not be merged on the second genome
  numberOfGenesInsideOfNonMergedIntervalsG2 = 0

  # The script has two steps:
  # 1st step: search for intervals that are present on the first species
  # 2nd step: search for intervals that are present on the second species

  # First Step ################################################################

  # Register the begin of the first step on the log file
  if (!is.null(logFileName)) {
    cat("Step 1 - Searching intervals on the first species\n", file=logFileName)
  }

  # Get the list of chromosomes of the first species
  chromosomeList = unique(as.character(orthologyTable$c1))

  # Process the list of genes of each chromosome
  for (chr in chromosomeList) {

    # Get just the genes of the chromosome chr
    tab = orthologyTable[orthologyTable$c1==chr,]

    # Create a table with the values of the begin and the end of the genes
    a = data.frame(inf=tab$inf1, sup=tab$sup1)

    # Make the union of the two sets of segments
    interval1 = unionIntervals(a, F)
    if (nrow(interval1) > 0) {
      names(interval1) = c("inf1", "sup1")

      # Create a vector full of zeros and nrow(tab) positions
      group1 = rep(0, nrow(tab))

      # Fill the group1 vector
      # - group1[i] = 0, if the gene i does not overlap any other gene
      # - group1[i] = X, if the gene is inside of an interval (X is the line of
      #               the table interval1 where is located the overlapping
      #               interval.
      for(i in 1:nrow(interval1)) {
        li1 = which(tab$inf1 >= interval1$inf1[i] & tab$sup1 <= interval1$sup1[i])
        if (length(li1) > 1) {
          group1[li1] = i
        }
      }

      # Get the lines of the overlapping intervals
      lesGroupes = unique(group1)
      lesGroupes = lesGroupes[xor(lesGroupes,0)]

      # Number of overlapping intervals
      numberOfOverlappingIntervalsG1 = numberOfOverlappingIntervalsG1 + length(lesGroupes)

      # Number of genes inside of the overlapping intervals
      numberOfGenesInsideOfOverlappingIntervalsG1 = numberOfGenesInsideOfOverlappingIntervalsG1 + sum(group1>0)

      # Process each one of the overlapping intervals
      groupeOK = vector()
      for (i in 1:length(lesGroupes)) {

        # Get the interval
        gp = lesGroupes[i]
        # Get the genes which are inside of the interval
        a = tab[group1==gp,]
        # Sort the genes of the interval by inf1
        a = a[order(a$inf1),]

        if (length(unique(a$c2)) == 1) {
          # In this case, all the overlapping genes from the first species
          # are inside of the same chromosome on the second species

          # Building the correspondent interval on the second species
          inf2 = min(a$inf2)
          sup2 = max(a$sup2)
          a2 = orthologyTable[orthologyTable$c2 == a$c2[1] &
            orthologyTable$sup2 >= inf2 & orthologyTable$inf2 <= sup2,]

          # Get the list of ortholog genes of the second species which are inside of the interval
          intersec = intersect(as.character(a$g2), as.character(a2$g2))

          if (length(intersec) == length(a2$g2) & length(intersec) == length(a$g2)) {
            # In this case, the same genes of the second species that appears on a,
            # also appears on a2

            # We have to verify if they have the same order and orientation
            pe = rank(a$inf2)
            peS = ifelse(a$strand1 == a$strand2, pe, -pe)

            if (all(peS == (1:length(peS))) | all(peS == (-(length(peS):1))) ) {

              # In this case, we have same order and orientation in the left
              # extremity of the interval

              # Add the interval on the group of good intervals
              groupeOK = c(groupeOK, gp)

              # Add a merged gene
              nom1 = paste("m", a$g1[1], sep="")
              merge.g1 = c(merge.g1, nom1)
              merge.c1 = c(merge.c1, as.character(a$c1[1]))
              merge.inf1 = c(merge.inf1, interval1$inf1[gp])
              merge.sup1 = c(merge.sup1, interval1$sup1[gp])
              merge.strand1 = c(merge.strand1, 1)
              nom2 = paste("m", a$g2[1], sep="")
              merge.g2 = c(merge.g2, nom2)
              merge.c2 = c(merge.c2, as.character(a$c2[1]))
              merge.inf2 = c(merge.inf2, inf2)
              merge.sup2 = c(merge.sup2, sup2)
              merge.strand2 = c(merge.strand2, a$strand1[1]*a$strand2[1])

              # Write on the log file the names of the merged gene and of the
              # genes that are inside of the interval
              if (!is.null(logFileName)) {
                cat(nom1, length(a$g1), as.character(a$g1), sep="\t", file=logFileName, append=T)
                cat("\n", file=logFileName, append=T)
                cat(nom2, length(a$g2), as.character(a$g2), sep="\t", file=logFileName, append=T)
                cat("\n", file=logFileName, append=T)
              }

              # Increase the number of merged intervals
              numberOfMergedIntervalsG1 = numberOfMergedIntervalsG1 + 1

              # Increase the number of genes which are inside of merged intervals
              numberOfGenesInsideOfMergedIntervalsG1 = numberOfGenesInsideOfMergedIntervalsG1 + length(a$g1)

              # Remove the merged genes from the list of ortholog genes
              orthologyTable = orthologyTable[!is.element(orthologyTable$g1, a$g1),]
              genesOnMergedIntervals = c(genesOnMergedIntervals, as.character(a$g1))

            } else {

              # In this case, we try by looking to the right extremity
              # of the interval

              # Sort the genes of the interval by sup1
              a = a[order(a$sup1),]

              # We have to verify if they have the same order and orientation
              pe = rank(a$sup2)
              peS = ifelse(a$strand1 == a$strand2, pe, -pe)

              if (all(peS == (1:length(peS))) | all( peS ==(-(length(peS):1))) ) {

                # In this case, we have same order and orientation in the left
                # extremity of the interval

                # Add the interval on the group of good intervals
                groupeOK=c(groupeOK,gp)

                # Add a merged gene
                nom1 = paste("m", a$g1[1], sep="")
                merge.g1 = c(merge.g1, nom1)
                merge.c1 = c(merge.c1, as.character(a$c1[1]))
                merge.inf1 = c(merge.inf1, interval1$inf1[gp])
                merge.sup1 = c(merge.sup1, interval1$sup1[gp])
                merge.strand1 = c(merge.strand1, 1)
                nom2 = paste("m", a$g2[1], sep="")
                merge.g2 = c(merge.g2, nom2)
                merge.c2 = c(merge.c2, as.character(a$c2[1]))
                merge.inf2 = c(merge.inf2, inf2)
                merge.sup2 = c(merge.sup2, sup2)
                merge.strand2 = c(merge.strand2, a$strand1[1]*a$strand2[1])

                # Write on the log file the names of the merged gene and of the
                # genes that are inside of the interval
                if (!is.null(logFileName)) {
                  cat(nom1, length(a$g1), as.character(a$g1), sep="\t", file=logFileName, append=T)
                  cat("\n", file=logFileName, append=T)
                  cat(nom2, length(a$g2), as.character(a$g2), sep="\t", file=logFileName, append=T)
                  cat("\n", file=logFileName, append=T)
                }

                # Increase the number of merged intervals
                numberOfMergedIntervalsG1 = numberOfMergedIntervalsG1 + 1

                # Increase the number of genes which are inside of merged intervals
                numberOfGenesInsideOfMergedIntervalsG1 = numberOfGenesInsideOfMergedIntervalsG1 + length(a$g1)

                # Remove the merged genes from the list of ortholog genes
                orthologyTable = orthologyTable[!is.element(orthologyTable$g1, a$g1),]
                genesOnMergedIntervals = c(genesOnMergedIntervals, as.character(a$g1))
              }
            }
          }
        }
      }
      # Update the list of genes that are inside of intervals
      genesOnIntervals = c(genesOnIntervals, as.character(tab$g1[group1>0]))
    }
  }

  # Update the number of non merged intervals
  numberOfNonMergedIntervalsG1 = numberOfOverlappingIntervalsG1 - numberOfMergedIntervalsG1
  # Update the number of genes which are inside of non merged intervals
  numberOfGenesInsideOfNonMergedIntervalsG1 = numberOfGenesInsideOfOverlappingIntervalsG1 - numberOfGenesInsideOfMergedIntervalsG1

  # Second Step ###############################################################

  # Register the begin of the second step on the log file
  if (!is.null(logFileName)) {
    cat("Step 2 - Searching intervals on the second species\n", file=logFileName, append=T)
  }

  # Get the list of chromosomes of the second species
  chromosomeList = unique(as.character(orthologyTable$c2))

  # Process the list of genes of each chromosome
  for (chr in chromosomeList) {

    # Get just the genes of the chromosome chr
    tab = orthologyTable[orthologyTable$c2==chr,]

    # Create a table with the values of the begin and the end of the genes
    a = data.frame(inf=tab$inf2, sup=tab$sup2)

    # Make the union of the two sets of segments
    interval2 = unionIntervals(a, F)
    if (nrow(interval2) > 0) {
      names(interval2) = c("inf2", "sup2")

      # Create a vector full of zeros and nrow(tab) positions
      group2 = rep(0, nrow(tab))

      # Fill the group2 vector
      # - group2[i] = 0, if the gene i does not overlap any other gene
      # - group2[i] = X, if the gene is inside of an interval (X is the line of
      #               the table interval1 where is located the overlapping
      #               interval.
      for(i in 1:nrow(interval2)) {
        li2 = which(tab$inf2 >= interval2$inf2[i] & tab$sup2 <= interval2$sup2[i])
        if (length(li2) > 1) {
          group2[li2] = i
        }
      }

      # Get the lines of the overlapping intervals
      lesGroupes = unique(group2)
      lesGroupes = lesGroupes[xor(lesGroupes,0)]

      # Number of overlapping intervals
      numberOfOverlappingIntervalsG2 = numberOfOverlappingIntervalsG2 + length(lesGroupes)

      # Number of genes inside of the overlapping intervals
      numberOfGenesInsideOfOverlappingIntervalsG2 = numberOfGenesInsideOfOverlappingIntervalsG2 + sum(group2>0)

      # Process each one of the overlapping intervals
      groupeOK = vector()
      for (i in 1:length(lesGroupes)) {

        # Get the interval
        gp = lesGroupes[i]
        # Get the genes which are inside of the interval
        a = tab[group2==gp,]
        # Sort the genes of the interval by inf1
        a = a[order(a$inf2),]

        if (length(unique(a$c2)) == 1) {
          # In this case, all the overlapping genes from the first species
          # are inside of the same chromosome on the second species

          # Building the correspondent interval on the first species
          inf1 = min(a$inf1)
          sup1 = max(a$sup1)
          a1 = orthologyTable[orthologyTable$c1 == a$c1[1] &
            orthologyTable$sup1 >= inf1 & orthologyTable$inf1 <= sup1,]

          # Get the list of ortholog genes of the first species which are inside of the interval
          intersec = intersect(as.character(a$g2), as.character(a1$g2))

          if (length(intersec) == length(a1$g2) & length(intersec) == length(a$g2)) {
            # In this case, the same genes of the first species that appears on a,
            # also appears on a2

            # We have to verify if they have the same order and orientation
            pe = rank(a$inf1)
            peS = ifelse(a$strand1 == a$strand2, pe, -pe)

            if (all(peS == (1:length(peS))) | all(peS == (-(length(peS):1))) ) {

              # In this case, we have same order and orientation in the left
              # extremity of the interval

              # Add the interval on the group of good intervals
              groupeOK = c(groupeOK, gp)

              # Add a merged gene
              nom1 = paste("m", a$g1[1], sep="")
              merge.g1 = c(merge.g1, nom1)
              merge.c1 = c(merge.c1, as.character(a$c1[1]))
              merge.inf1 = c(merge.inf1, inf1)
              merge.sup1 = c(merge.sup1, sup1)
              merge.strand1 = c(merge.strand1, 1)
              nom2 = paste("m", a$g2[1], sep="")
              merge.g2 = c(merge.g2, nom2)
              merge.c2 = c(merge.c2, as.character(a$c2[1]))
              merge.inf2 = c(merge.inf2, interval2$inf2[gp])
              merge.sup2 = c(merge.sup2, interval2$sup2[gp])
              merge.strand2 = c(merge.strand2, a$strand1[1]*a$strand2[1])

              # Write on the log file the names of the merged gene and of the
              # genes that are inside of the interval
              if (!is.null(logFileName)) {
                cat(nom1, length(a$g1), as.character(a$g1), sep="\t", file=logFileName, append=T)
                cat("\n", file=logFileName, append=T)
                cat(nom2, length(a$g2), as.character(a$g2), sep="\t", file=logFileName, append=T)
                cat("\n", file=logFileName, append=T)
              }

              # Increase the number of merged intervals
              numberOfMergedIntervalsG2 = numberOfMergedIntervalsG2 + 1

              # Increase the number of genes which are inside of merged intervals
              numberOfGenesInsideOfMergedIntervalsG2 = numberOfGenesInsideOfMergedIntervalsG2 + length(a$g1)

              # Remove the merged genes from the list of ortholog genes
              orthologyTable = orthologyTable[!is.element(orthologyTable$g1, a$g1),]
              genesOnMergedIntervals = c(genesOnMergedIntervals, as.character(a$g1))

            } else {

              # In this case, we try by looking to the right extremity
              # of the interval

              # Sort the genes of the interval by sup2
              a = a[order(a$sup2),]

              # We have to verify if they have the same order and orientation
              pe = rank(a$sup1)
              peS = ifelse(a$strand1 == a$strand2, pe, -pe)

              if (all(peS == (1:length(peS))) | all( peS ==(-(length(peS):1))) ) {

                # In this case, we have same order and orientation in the left
                # extremity of the interval

                # Add the interval on the group of good intervals
                groupeOK=c(groupeOK,gp)

                # Add a merged gene
                nom1 = paste("m", a$g1[1], sep="")
                merge.g1 = c(merge.g1, nom1)
                merge.c1 = c(merge.c1, as.character(a$c1[1]))
                merge.inf1 = c(merge.inf1, inf1)
                merge.sup1 = c(merge.sup1, sup1)
                merge.strand1 = c(merge.strand1, 1)
                nom2 = paste("m", a$g2[1], sep="")
                merge.g2 = c(merge.g2, nom2)
                merge.c2 = c(merge.c2, as.character(a$c2[1]))
                merge.inf2 = c(merge.inf2, interval2$inf2[gp])
                merge.sup2 = c(merge.sup2, interval2$sup2[gp])
                merge.strand2 = c(merge.strand2, a$strand1[1]*a$strand2[1])

                # Write on the log file the names of the merged gene and of the
                # genes that are inside of the interval
                if (!is.null(logFileName)) {
                  cat(nom1, length(a$g1), as.character(a$g1), sep="\t", file=logFileName, append=T)
                  cat("\n", file=logFileName, append=T)
                  cat(nom2, length(a$g2), as.character(a$g2), sep="\t", file=logFileName, append=T)
                  cat("\n", file=logFileName, append=T)
                }

                # Increase the number of merged intervals
                numberOfMergedIntervalsG2 = numberOfMergedIntervalsG2 + 1

                # Increase the number of genes which are inside of merged intervals
                numberOfGenesInsideOfMergedIntervalsG2 = numberOfGenesInsideOfMergedIntervalsG2 + length(a$g1)

                # Remove the merged genes from the list of ortholog genes
                orthologyTable = orthologyTable[!is.element(orthologyTable$g1, a$g1),]
                genesOnMergedIntervals = c(genesOnMergedIntervals, as.character(a$g1))
              }
            }
          }
        }
      }
      # Update the list of genes that are inside of intervals
      genesOnIntervals = c(genesOnIntervals, as.character(tab$g1[group2>0]))
    }
  }

  genesOnIntervals = unique(genesOnIntervals)

  numberOfNonMergedIntervalsG2 = numberOfOverlappingIntervalsG2 - numberOfMergedIntervalsG2
  numberOfGenesInsideOfNonMergedIntervalsG2 = numberOfGenesInsideOfOverlappingIntervalsG2 - numberOfGenesInsideOfMergedIntervalsG2

  # Remove the genes that were inside of intervals but could not be merged
  orthologyTable = orthologyTable[!is.element(orthologyTable$g1, genesOnIntervals),]

  # Create a table with the merged intervals
  mergedIntervals = data.frame(g1 = merge.g1,
                               c1 = merge.c1,
                               inf1 = merge.inf1,
                               sup1 = merge.sup1,
                               strand1 = merge.strand1,
                               g2 = merge.g2,
                               c2 = merge.c2,
                               inf2 = merge.inf2,
                               sup2 = merge.sup2,
                               strand2 = merge.strand2)

  # Combine the table of non overlapping genes with the table of
  # merged intervals
  nonOverlappingGenesAndMergedIntervals = rbind(orthologyTable, mergedIntervals)

  # Create a list with statistics about the intervals
  statistics = list(numberOfOverlappingIntervalsG1 = numberOfOverlappingIntervalsG1,
                    numberOfGenesInsideOfOverlappingIntervalsG1 = numberOfGenesInsideOfOverlappingIntervalsG1,
                    numberOfMergedIntervalsG1 = numberOfMergedIntervalsG1,
                    numberOfGenesInsideOfMergedIntervalsG1 = numberOfGenesInsideOfMergedIntervalsG1,
                    numberOfNonMergedIntervalsG1 = numberOfNonMergedIntervalsG1,
                    numberOfGenesInsideOfNonMergedIntervalsG1 = numberOfGenesInsideOfNonMergedIntervalsG1,
                    numberOfOverlappingIntervalsG2 = numberOfOverlappingIntervalsG2,
                    numberOfGenesInsideOfOverlappingIntervalsG2 = numberOfGenesInsideOfOverlappingIntervalsG2,
                    numberOfMergedIntervalsG2 = numberOfMergedIntervalsG2,
                    numberOfGenesInsideOfMergedIntervalsG2 = numberOfGenesInsideOfMergedIntervalsG2,
                    numberOfNonMergedIntervalsG2 = numberOfNonMergedIntervalsG2,
                    numberOfGenesInsideOfNonMergedIntervalsG2 = numberOfGenesInsideOfNonMergedIntervalsG2)

  # Create a table of genes that were discarded (overlapping genes that
  # are inside of intervals that could not be merged)
  discardedGenes = genesOnIntervals[!is.element(genesOnIntervals, genesOnMergedIntervals)]

  # Create the list of genes (one list for each type)
  genesList = list(genesOnIntervals = genesOnIntervals,
                   genesOnMergedIntervals = genesOnMergedIntervals,
                   discardedGenes = discardedGenes)

  # Return the processed data
  # 1 - Table that contains the genes that were not inside of intervals
  #     and the merged intervals
  # 2 - Table that contains the genes that were not inside of intervals
  # 3 - Table that contains the merged intervals
  # 4 - Vector with statistics values
  # 5 - List of genes
  return (list(nonOverlappingGenesAndMergedIntervals = nonOverlappingGenesAndMergedIntervals,
               nonOverlappingGenes = orthologyTable,
               mergedIntervals = mergedIntervals,
               statistics = statistics,
               genesList = genesList))
}
