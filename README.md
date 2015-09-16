# StrainTypeMer

StrainTypeMer is a tool for comparing strains of bacteria using their kmer content.  StrainTyper is multithreaded
and can perform comparisons within minutes.  Fast-mode compares strains in seconds by comparing a subset of the kmers.

The results of StrainTypeMer include:
* Simlarity Table - displays the percentage of kmers shared between each strain combination
* Distance Matrix - png that used clustering to show the strain relationships
* Sample Summary - written to standard
    1. Estimated genome size of each strain
    2. Coverage cutoff
    3. ARDB - Antibiotic Resistance Database content
        * Each strain is queried for kmers for gene present in the ARDB
        * Number of kmers are report for each gene family


<hr>

### CONSIDERATIONS

In order to use StrainTypeMer you must first count kmers using [jellyfish 2.0 of greater](https://github.com/gmarcais/Jellyfish)
for each strain.  While preprocessing is not required it can be performed prior to kmer counting.  We generally remove
adapters and low quality/reads prior to kmer counting [SHOW DATA BEFORE / AFTER].

#### Kmer counting
The kmer size can be varied, we used 31 bp kmers to validate StrainTypeMer.
It is important to keep the hash size equal when count each strain.


`jellyfish count -m 31 -L 2 -t 6 -o strain.jf -s 4G -C strain_1.fq strain_2.fq `

the -L flag prevents jellyfish outputing kmers seen only once, it is not required. However these kmers will be excluded
anyway, the file size will be reduced and the StrainTypeMer will not iterate over the kmer.

