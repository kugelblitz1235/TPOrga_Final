#!/bin/bash
for f in genomes/*_RefSeq.fasta; do
    name=$(basename "$f" .fasta)
    mkdir -p "reads/$name"
    ./mason illumina -n 100 -N 10 -f $f -o "reads/${name}/reads.fasta"
    split -l 2 -d -a 2 --additional-suffix=".fasta" "reads/${name}/reads.fasta" "reads/${name}/"
    rm "reads/${name}/reads.fasta" "reads/${name}/reads.fasta.sam"
done
