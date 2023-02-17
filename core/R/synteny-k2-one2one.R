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
# This file contains the definition of the functions filterOrthologyTableK2 and
# buildK2Graph.
#
# Author: Claire Lemaitre
# Reviewed by: Christian Baudet
# #############################################################################

# #############################################################################
# FUNCTION filterOrthologyTableK2(orthologyTable)
#
# About this function:
# 1 - This function receives a table of ortholog genes which has the columns:
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
# 2 - The columns of the table must be named and ordered as described above.
#
# 3 - The input table can not have duplication of genes and must have
#     only one to one relationships between the genes of the two species.
#
# 4 - The funtion will process this table and return just the pairs of ortholog
#     genes which do not show conflict of type I or II. It also removes all
#     genes that do not have a correspondent ortholog gene.
#
# 5 - The function will return a dataframe that has the same fields of the
#     received dataframe.
#
# #############################################################################
filterOrthologyTableK2 = function(orthologyTable) {

  # Change the column names
  names(orthologyTable) = c("g1", "c1", "inf1", "sup1", "strand1", "g2", "c2", "inf2", "sup2", "strand2")

  # Build the graph
  graph = buildK2Graph(orthologyTable)

  # Build a list of unique gene names
  gene1Names = unique(c(as.character(graph$anchorAgene1), as.character(graph$anchorBgene1)))

  # Filter the orthology table to preserve just the pairs (g1, g2) if g1 is
  # on the list of names that appears on the graph
  filteredOrthologyTable = orthologyTable[is.element(as.character(orthologyTable$g1), gene1Names),]

  # Sort the filtered orthology table by the chromosome and start position of
  # the gene on the first species
  filteredOrthologyTable = filteredOrthologyTable[order(filteredOrthologyTable$c1, filteredOrthologyTable$inf1),]

  # Create a vector that says, for each line of the filtered orthology table,
  # the rank of the row when we sort the table by the chromosome and the
  # begining of the gene on the second species
  per = order(order(filteredOrthologyTable$c2, filteredOrthologyTable$inf2))

  # Creates a vector that says, for each line of the filtered orthology table,
  # the signed rank of the row when we sort the table by the chromosome and the
  # begining of the gene on the second species. The rank value will be
  # positive if the genes of the two species have the same direction and
  # negative, otherwise
  perS = ifelse(filteredOrthologyTable$strand1==filteredOrthologyTable$strand2, per, -per)

  # Create a vector that has the result of:
  # distance[i] = perS[i+1] - perS[i]
  # It indicates the distance between two consecutive rows
  distance = diff(perS)

  # Create a vector which has the number of the lines on the vector
  # distance that have distance 1
  distance1 = which(distance==1)

  # Create a new dataframe that have the relationships between
  # every row that has distance 1
  auxDataFrame1 = data.frame(filteredOrthologyTable[distance1, c("g1", "c1", "c2")],
                             filteredOrthologyTable[distance1+1, c("g1", "c1", "c2")])

  # Make the chromosome matching
  auxDataFrame2 = auxDataFrame1[auxDataFrame1$c1==auxDataFrame1$c1.1 & auxDataFrame1$c2==auxDataFrame1$c2.1,]

  # Build new a list of unique gene names
  gene1Names = unique(c(as.character(auxDataFrame2$g1), as.character(auxDataFrame2$g1.1)))

  # Build the final result. A data frame that has only the genes that
  # appears on the list gene1Names
  toReturn = filteredOrthologyTable[is.element(as.character(filteredOrthologyTable$g1), gene1Names),]

  return (toReturn)
}

# #############################################################################
# FUNCTION buildK2Graph(orthologyTable)
#
# About this function:
# 1 - This function receives a table of ortholog genes which has the columns:
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
# 2 - The columns of the table must be named and ordered as described above.
#
# 3 - The input table can not have duplication of genes and must have
#     only one to one relationships between the genes of the two species.
#
# 4 - This table is used to define a set of anchors. An anchor can be define
#     as an object (c1, c2, r1, r2, s), where:
#     c1 - Chromosome where the gene g1 is located;
#     c2 - Chromosome where the gene g2 is located;
#     r1 - Rank of the gene g1 in the chromosome c1;
#     r2 - Rank of the gene g2 in the chromosome c2;
#     s  - s = 1 if strand1==strand2 and s = -1, otherwise.
#     (Note that rank is the number of the gene on the chromosome).
#
# 5 - This function will build a graph of anchors. An arc links two anchors
#     A and B if they are colinear and the distance between them is smaller
#     than a parameter k. (Note that the distance is measured with their
#     rank positions, and not by their position on the chromosome)
#
# 6 - Un arc AB will be present on the graph if:
#     - A.r1 < B.r1
#     - d(A,B) <= k
#     - (A.r2 < B.r2 and A.s = B.s = 1) or (A.r2 > B.r2 and A.s = B.s = -1)
#
# 7 - This function uses k = 2
#
# 8 - The function will return a data frame that contains the columns:
#     anchorAgene1 - Name of the gene on the first species that is related
#                    with the anchor A of an arc AB which is present into
#                    the graph
#     anchorBgene1 - Name of the gene on the first species that is related
#                    with the anchor B of an arc AB which is present into
#                    the graph
#     label        - label that identifies the type of the arc:
#                    ok      - consecutive genes
#                    inser1  - Insertion of a gene on the first genome
#                    inser2  - Insertion of a gene on the second genome
#                    inser12 - Insertion of genes in both genomes
#
# #############################################################################
buildK2Graph = function(orthologyTable) {

  # Auxiliary variables
  consecutiveGenes=NULL
  insertionOnG2=NULL
  insertionOnG1=NULL
  insertionOnBoth=NULL

  # Create a table sorted by the name of the chromosome and the
  # begining position of the gene (on the first species)
  tab1 = orthologyTable[order(orthologyTable$c1, orthologyTable$inf1),]

  # Create a table sorted by the name of the chromosome and the
  # begining position of the gene (on the second species)
  tab2 = orthologyTable[order(orthologyTable$c2, orthologyTable$inf2),]

  # Create a vector that says, for each line of the table tab1, the rank of
  # the row when we sort the table by the chromosome and the begining of
  # the gene on the second species
  per = order(order(tab1$c2, tab1$inf2))

  # Create a vector that says, for each line of the table tab1, the signed
  # rank of the row when we sort the table by the chromosome and the
  # begining of the gene on the second species. The rank value will be
  # positive if the genes of the two species have the same direction and
  # negative, otherwise
  perS = ifelse(tab1$strand1==tab1$strand2, per, -per)

  # Create a vector that has the result of:
  # distance1[i] = perS[i+1] - perS[i]
  # It indicates the distance between two consecutive rows
  distance1 = diff(perS)

  # Create a vector which has the number of the lines on the vector
  # distance1 that have distance 1
  distance1k1= which(distance1 == 1)

  # Create a vector which has the number of the lines on the vector
  # distance1 that have distance 2
  distance1k2 = which(distance1 == 2)

  if (length(distance1k1) > 0) {
    # Consecutive genes
    # Create a data frame that has the consecutive genes (distance = 1) (ok label)
    consecutiveGenes = data.frame(tab1[distance1k1, c("g1", "c1", "c2")], tab1[distance1k1+1, c("g1", "c1", "c2")], lab="ok")
  }

  if (length(distance1k2) > 0) {
    # In this case, we have a gene inserted on the genome of the second species
    # Create a data frame with the genes that have distance 2 (inser2 label)
    insertionOnG2 = data.frame(tab1[distance1k2, c("g1", "c1", "c2")], tab1[distance1k2+1, c("g1", "c1", "c2")], lab="inser2")
  }

  # Create a vector that has the result of:
  # distance2[i] = perS[i+2] - perS[i]
  distance2 = diff(perS, lag=2)

  # Create a vector with the lines of the vector distance2 that have
  # distance2 = 1 (excluding the adjacent genes)
  distance2k1 = which(distance2==1 & distance1[-length(distance1)] != 1)

  if (length(distance2k1) > 0) {
    # In this case, we have a gene inserted on the genome of the first species
    # Create a data frame with the genes (inser1 label)
    insertionOnG1 = data.frame(tab1[distance2k1, c("g1", "c1", "c2")], tab1[distance2k1+2, c("g1", "c1", "c2")], lab="inser1")
  }

  # Create a vector with the lines of the vector distance2 that have
  # distance2 = 2 (excluding the adjacent genes)
  distance2k2 = which(distance2==2 & distance1[-length(distance1)] != 1)

  if (length(distance2k2) > 0) {
    # In this case, we have a gene inserted on both genomes
    # Create a data frame with the genes (inser12 label)
    insertionOnBoth = data.frame(tab1[distance2k2, c("g1", "c1", "c2")], tab1[distance2k2+2, c("g1", "c1", "c2")], lab="inser12")
  }

  # Bind all produced data frames (at this moment we are not looking if
  # the chromosomes match)
  notFilteredGeneGraph = rbind(consecutiveGenes, insertionOnG2, insertionOnG1, insertionOnBoth)

  # Perform the chromosome matching
  geneGraph = notFilteredGeneGraph[notFilteredGeneGraph$c1==notFilteredGeneGraph$c1.1 & notFilteredGeneGraph$c2==notFilteredGeneGraph$c2.1,]

  # Rename the columns
  names(geneGraph)[names(geneGraph)=="g1"] = "anchorAgene1"
  names(geneGraph)[names(geneGraph)=="g1.1"] = "anchorBgene1"
  names(geneGraph)[names(geneGraph)=="lab"] = "label"

  # Return the final data frame
  return (geneGraph)
}
