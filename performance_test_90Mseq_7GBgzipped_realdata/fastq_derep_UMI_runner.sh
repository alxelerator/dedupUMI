#!/bin/bash

# runner for UMI dereplicator DEMO
# a.bossers@uu.nl // alex.bossers@wur.nl

date

./deduplicate_UMI.pl \
		--input-fastq1 test_R1.fastq.gz \
		--input-fastq2 test_R2.fastq.gz \
		--input-fastq3 test_R3.fastq.gz \
		--output-fastq1 derep_fullR1.fq.gzip \
		--output-fastq2 derep_fullR2.fq.gzip \
		--output-fastq3 derep_fullR3.fq.gzip \

date
