# StrainTypeMer
<hr>
##### StrainTypeMer is a tool for pairwise comparison multiple genomes. Comparisons are made be comparing the 31bp kmers found between different strains/samples.

##### Features include:
* Reference free comparison
* pan-genome and core-genome comparison using genome filter (see below)
* Pairwise comparison of each input strain
* Filters poor quality sequence
* QC metrics
    * Estimated genome size
    * Kmer cutoff
    * Low coverage warning
* Detection of epidemiological important genes
* reporting of MLST type
* Clustering and dendrogram creation
* PDF output
* Text output of similarity table and QC metrics

<hr>
##### Inputs:
* input can be fastq, fasta, fastq.gz, or fasta.gz files as input.

<hr>
##### Installation:

1 If not already on your system install python 2.7 (I use anaconda)

2 StrainTypeMer requires Jellyfish installed with the python bindings.
``` bash
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz
tar xvf jellyfish-2.2.5.tar.gz
cd jellyfish-2.2.5
./configure --prefix=$HOME --enable-python-binding
make
make install

```
3 Install StrainTypeMer
```bash
git clone
cd StrainTypeMer
python setup.py install
```

<hr>
##### Usage:
StrainTypeMer has the following sub-commands:
 * update_mlst - This updates the MLST resources required to do MLST profiling
 * compare - Main command that take sequence files performs kmer filtering, MLST profiling, Gene detection, and comparison
    * Take >2 strains as input
 * count - a subset of compare that peforms kmer filtering, MLST profiling, gene detection, and writes the strain record to json format.
    * This command takes a single strain
    * TODO this will allow for client side pairwise comparison

##### Considerations:
StrainTypeMer was tested on Illumina datasets, where we found a minimum of 20X coverage was need for optimal strain
comparsion.



