#!/bin/bash

# runner for UMI dereplicator DEMO
# a.bossers@uu.nl // alex.bossers@wur.nl

# This example also shows that variable input output gzipped options can be used. They are handled based on given finput/output filenames.
echo -e "\n>>> Run test for R3-system (R1 r2 R3 input/output files having UMI in R2) <<<\n"
./deduplicate_UMI \
		--input-fastq1 ./demo/aR1.fq.gz \
		--input-fastq2 ./demo/aR2.fq \
		--input-fastq3 ./demo/aR3.fq \
		--output-fastq1 ./outR1.fq.gz \
		--output-fastq2 ./outR2.fq \
		--output-fastq3 ./outR3.fq 
	 	# add --verbose to output sequences to STDOUT


echo -e "\n>>> Run test for UMI in headers (new format) <<<\n"
#Test using UMI in the HEADERS instead of in R1 R2 R3 system
./deduplicate_UMI \
		--input-fastq1 ./demo/bR1.fq \
		--input-fastq2 ./demo/bR2.fq \
		--output-fastq1 ./out_headR1.fq \
		--output-fastq2 ./out_headR2.fq 

