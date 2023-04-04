# Cassis: Detection and refinment of genomic rearrangement breakpoints

Revival of Cassis, previous web page : [http://pbil.univ-lyon1.fr/software/Cassis/](http://pbil.univ-lyon1.fr/software/Cassis/)


## User manual

### Installation 

Requirements :

* Perl
* R
* Ghostscript
* [lastz](https://github.com/lastz/lastz)

Note: lastz binary must be put in directory [./lastz/](./lastz/)

### Usage

```
perl cassis.pl [options] <table.txt> <type> <dirGR> <dirGO> <outputdir>
```

or (new since version 2.0)

```
perl cassis.pl [options] <table.txt> <type> <GR.fasta> <GO.fasta> <outputdir>
```

#### Mandatory parameters

* `table` : tab separated file with orthologous genes or synteny blocks coordinates in both genomes (see expected format in [cassis_manual.pdf](cassis_manual.pdf) p.9)
* `type` : type of table, either `G` for genes or `B` for blocks.
* `dirGR` and `dirGO` : 1 directory per genome, each contains several fasta files (one for each chromosome or scaffold)
* `GR.fasta` and `GO.fasta` **[new in version 2.0]** : alternatively to dirs and numerous fasta files, the sequences of both genomes can be given as multi-fasta files
* `outputdir` : output directory (see the different output files in [cassis_manual.pdf](cassis_manual.pdf) p.11-15)

Detailed user manual of version 1.0, with description of all other optional parameters in [cassis_manual.pdf](cassis_manual.pdf)

### Novel features in version 2.0

* lastz level 1
* instead of giving the genome sequences as numerous unique-entry fasta files for each chromosome, we can now give one single multi-fasta file for each genome.




## Citation

- [Cassis: Detection of genomic rearrangement breakpoints.](http://bioinformatics.oxfordjournals.org/cgi/content/short/26/15/1897)
C. Baudet, C. Lemaitre, Z. Dias, C. Gautier, E. Tannier, M.-F. Sagot.
Bioinformatics, 2010 26(15):1897-1898.

- [Precise detection of rearrangement breakpoints in mammalian chromosomes.](http://www.biomedcentral.com/1471-2105/9/286) Lemaitre, C., Tannier, E., Gautier, C., and Sagot, M.-F. BMC Bioinformatics, 2008, 9(286), 15 pages. 