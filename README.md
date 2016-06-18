# StrainTypeMer
<hr>
##### StrainTypeMer is a tool for pairwise comparison multiple genomes. StrainTypeMer compares the 31bp kmers found between different strains/samples.

##### Features include:
* Reference free comparison
* coregenome or non-coregenome comparison using filter fasta file.
* Pairwise comparison of each input strain
* Filters poor quality sequence
* QC metrics
    * Estimated genome size
    * Kmer count cutoff
    * Low coverage warning
* Detection of epidemiological important genes [BETA `-ard` argument]
* MLST type
* Clustering and dendrogram creation
* PDF output (strain histograms and clustered similarity matrix)
* Text output of similarity table and QC metrics

<hr>
##### Inputs:
* input files can be fastq, fasta, fastq.gz, or fasta.gz

<hr>
##### Installation:

1.	If not already on your system install python 3.5
 * I suggest using anaconda for fewer headaches

2.	Install Jellyfish with the python bindings.
``` bash
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz
tar xvf jellyfish-2.2.5.tar.gz
cd jellyfish-2.2.5
./configure --prefix=$HOME --enable-python-binding
make
make install

```
3.	Install StrainTypeMer
```bash
git clone
cd StrainTypeMer
python setup.py install
```

<hr>
##### Usage:
full usage can be found by using:

```bash
python ~/StrainTypeMer/straintypemer/StrainTypeMer.py compare -h
```
StrainTypeMer has the following sub-commands:
 * update_mlst - This updates the MLST resources required to do MLST profiling
 * compare - Main command that take sequence files and performs kmer filtering, MLST profiling, Gene detection, and comparison
    * Take > 2 strains as input

Example command:

`python ~/StrainTypeMer/straintypemer/StrainTypeMer.py compare -t 3 --qual_score 15 -ard A10A.fq.gz:A10A A11A.fq.gz:A11A`


##### Considerations:
StrainTypeMer was tested on Illumina data sets, where we found a minimum of 20X coverage was need for optimal strain
comparison.

We are currently still developing StrainTypeMer and will continue to make improvements.





