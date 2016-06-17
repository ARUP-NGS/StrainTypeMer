#!/usr/bin/env bash

genome_dir="/Users/331-SimmonkLPTP/Documents/data/strain_typing/ancinetobacter_baumannii/genomes/"
count_dir="/Users/331-SimmonkLPTP/Documents/data/strain_typing/ancinetobacter_baumannii/counts_31/"
jellyfish="/Users/331-SimmonkLPTP/bin/jellyfish-2.2.5/bin/jellyfish"

for f in $genome_dir/*
    do echo $f;

    gzip -dc $f | $jellyfish count -m 31 -s "500M" -t 4 -C -o ${count_dir}$(basename $f ".fa.gz")_31.jf /dev/fd/0
done
