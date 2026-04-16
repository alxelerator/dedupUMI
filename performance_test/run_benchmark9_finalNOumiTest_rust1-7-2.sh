#!/bin/bash


# run benchmark on dedupe_UMI for RUST and perl version 1.4
# a.bossers@uu.nl
# input file roughly 70M sequences R2 system gzipped


echo -e "\n\nTesting RUST172 implementation NO UMI use test...\n"

./deduplicate_UMI.rust172 --version
echo

./deduplicate_UMI.rust172 \
   --input-fastq1 NGS-26-CA-s10_rust_R1.fastq.gz \
   --input-fastq2 NGS-26-CA-s10_rust_R2.fastq.gz \
   --output-fastq1 NGS-26-CA-s10_dedupeRust_R1.fq.gz \
   --output-fastq2 NGS-26-CA-s10_dedupeRust_R2.fq.gz \
   --noUMI \
   --timing

date
echo "All done!"
