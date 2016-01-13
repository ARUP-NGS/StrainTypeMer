# StrainTypeMer

StrainTypeMer is a tool for comparing bacterial strains. StrainTypeMer uses kmer count output from jellyfish to 
compare strains. StrainTyper is multithreaded and can perform comparisons within minutes.  Fast-mode compares 100,000 
kmers and produces results nearly identical to a full comparison in less than a minute.

The results of StrainTypeMer include:

* Simlarity Table - displays the percentage of kmers shared between each strain combination
* Distance Matrix - png that used clustering to show the strain relationships
* Histogram Plots - PDF showing results of five strains per page.
* Sample results - written to standard out
    1. Estimated genome size of each strain
    2. Coverage cutoff
    3. MLST profile (for Acinetobacter;
    4. ARDB - Antibiotic Resistance Database content
        * Each strain is queried for kmers for gene present in the ARDB
        * Number of kmers are report for each gene family


<hr>
### REQUIREMENTS

In order to use StrainTypeMer you must first count kmers using [jellyfish 2.2.3 of greater](https://github.com/gmarcais/Jellyfish)
for each strain.  While preprocessing is not required it can be performed prior to kmer counting.  We generally remove
adapters and low quality/reads prior to kmer counting. [SHOW DATA BEFORE / AFTER]. 

Future version may take fastq files as input.

#### Kmer counting
The kmer size can be varied, we used 31 bp kmers to validate StrainTypeMer.
It is important to keep the hash size equal when count each strain. 


`jellyfish count -m 31 -L 2 -t 6 -o strain.jf -s 4G -C strain_1.fq strain_2.fq `

the -L flag prevents jellyfish outputing kmers seen only once, it is not required. However these kmers will be excluded
anyway, the file size will be reduced and the StrainTypeMer will not iterate over the kmer.

### DEPENDENCIES

* StrainTyper requires that `Jellyfish 2.2.3` be installed with the python swig bindings enabled. 
* Python 2.7  <em> - The anaconda distribution of python includes Matplotlib and Numpy
    * Matplotlib
    * Numpy





