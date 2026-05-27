# Deduplicate FASTQ sequences using UMIs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20412197.svg)](https://doi.org/10.5281/zenodo.20412197)

Dedupuplicate_UMI is a tool written in *rust* to remove **alignment-free** exact-duplicate FASTQ read-pairs using the UMI sequences-library approach. 
Duplicates are identified using a concatenation of R1, R2 and UMI sequences approach. For duplicate 
molecules the read-pair with the highest **total base PHRED quality score** is retained.

## Why deduplicate?

Exact duplicate reads (including identical UMI sequences) are typically the result of over-amplification during PCR in the sequencing library preparation. For quantification 
and some assembly-based approaches you want these overamplified sequences removed.  

However, identical read sequences without the same UMI can occur naturally and may represent independent molecules that were sequenced multiple 
times. These should generally **NOT** be removed. The probability of observing such biologically relevant duplicates increases with deeper sequencing. Typically the occurrence 
of erroneous sequence duplicates is linked to input DNA concentrations and the number of amplification cycles used in the library prep.  

Exact duplicates should only be removed when both the read sequences **and** the UMI sequence are identical (amplification artefact).

## Duplicate definition

Reads are considered duplicates when the following combination is identical:
 `R1 sequence + R2 sequence + UMI`

*For the older R3 system this becomes: R1 sequence + R2 sequence + R3 sequence*


## Input and output

NEW system R1 R2 **UMI in header**:
- Input: R1 and R2 FASTQ readset. UMI should be present in the header of at least R1!
- Output: R1 and R2 filtered for EXACT duplicates based on concatenated seq of R1 R2 keeping highest TOTAL qual score.

Older system R1 R2 R3 **UMI in R2**:
- Input: R1 R2 and R3 FASTQ readset. For the best total quality to work the UMI should be present in R2!
- Output: R1 R2 R3 filtered for EXACT duplicates based on concatenated seq of R1 R2 R3 keeping highest TOTAL qual score.


## Important notes

- Binary rust build based on `x86_64-unknown-linux-musl` for maximum HPC compatibility.
- It writes out the last sequences found of a duplicate set having the highest TOTAL qualityscore.
- Input FASTQ files must follow the standard 4-line FASTQ format, starting at the first record.
  + Sequences in R1 R2 (and R3) should be in same order and NOT INTERLEAVED!!
- UMI location (header or separate FASTQ file) is detected automatically unless `--noUMI` is used
- If using `--noUMI` we just filter exact duplicates. This most likely filters too harsh since 
  exact duplicates can occur in natural good quality datasets. Especially at high sequencing depths.
- Gzip I/O is now handled natively via `rust::flate2` (zlib-ng backend).
  + No external gzip/zcat/pigz dependency required anymore from version 1.4r and up.
  + Gzip is the only parallelised part of the tool and its greedy for the number of cores. Limit the number of cores using `--max-cores` when needed.


## Requirements

- FASTQ files in fixed 4-line format having read-pairs in separate files (DO NOT use INTERLEAVED FASTQ files)!
- Since the rust version 1.4r Gzip I/O is now handled natively via rust::flate2 (zlib-ng backend).
  No external gzip/zcat/pigz dependency required anymore.


## Command line usage

Run `./deduplicateUMI --help` to see the built-in help.  
Do not forget to give deduplicate_UMI execution rights (for convenience): `chmod u+x deduplicate_UMI`

```
DedupUMI version 1.7.3 (2026-04-16)
https://github.com/alxelerator/dedupUMI
Binary build for HPC compatibility against x86_64-unknown-linux-musl
a.bossers@uu.nl / alex.bossers@wur.nl

Usage:
  deduplicate_UMI \
      --input-fastq1 <R1.fastq[.gz]> \
      --input-fastq2 <R2.fastq[.gz]> \
      --output-fastq1 <R1_out.fastq[.gz]> \
      --output-fastq2 <R2_out.fastq[.gz]> \
      [--input-fastq3 <UMI.fastq[.gz]>] \
      [--output-fastq3 <UMI_out.fastq[.gz]>] \
      [--output_counts <counts.tab>] \
      [--noUMI] \
      [--max-cores] \
      [--timing] \
      [--verbose] 

Required arguments:
  --input-fastq1 <file>     FASTQ read1 input file (plain or gz)
  --input-fastq2 <file>     FASTQ read2 input file (plain or gz)
  --output-fastq1 <file>    FASTQ read1 output file (plain or gz)
  --output-fastq2 <file>    FASTQ read2 output file (plain or gz)

Optional arguments:
  --input-fastq3 <file>     UMI FASTQ file (R3 system). If omitted, UMI is
                            extracted from the read header
  --output-fastq3 <file>    Output FASTQ for R3 system UMIs (required if input-fastq3 given)
  --output_counts <file>    Append input/output read counts to output tabular file
  --noUMI                   Two-file mode without UMI. Deduplicate on R1+R2 only (see README)
  --max-cores               By default gzip compression is greedy using all available cores.
  --timing                  Print timing information to monitor the two major steps: read-hash and gz write.
  --verbose                 Print sequences to STDOUT (debug output)
  --version                 Print version information
  --help                    Show this help message

Notes:
  - UMI location (header or separate FASTQ file) is detected automatically unless --noUMI is used
  - If using --noUMI we just filter exact duplicates. This most likely filters too harsh since
    exact duplicates can occur in natural good quality datasets. Especially at high sequencing depths.
  - Gzip I/O is now handled natively via rust::flate2 (zlib-ng backend).
    No external gzip/zcat/pigz dependency required anymore from version 1.4r and up.
    Gzip is the only parallelised part and its greedy for the number of cores. Limit the number of cores using --max-cores.

Examples:
  UMI in header:
    deduplicate_UMI --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz \
                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz
  R3 UMI system:
    deduplicate_UMI --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz --input-fastq3 s1_R3.fastq.gz \
                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz --output-fastq3 s1_R3.dedup.fastq.gz
  No UMI:
    deduplicate_UMI --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz --noUMI \
                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz
```



## How it works

DedupUMI was initially implemented as a lightweight Perl script using an in-memory
hash table keyed by concatenation of:
`R1_sequence + R2_sequence + UMI`
This perl version has been replaced by a RUST rewrite from version 1.4r and up!  

This allows duplicate detection without alignment, sorting, or external tools.
Because the algorithm only performs sequential FASTQ reading combined
with hash lookups, runtime is typically limited by disk I/O
(decompression, compression and writing FASTQ files), not by the
deduplication logic itself.

For each read-pair the script:
1. Extracts the UMI (from the FASTQ header or the R3 file)
2. Constructs a (hash)key: R1_sequence + R2_sequence + UMI
3. If the key is new → the sequence reads and its quality values are stored
4. If the key already exists → the stored sequence reads and its quality scores are replaced only if the
   new readset has a higher total base quality score

After all reads are processed the remaining unique reads are written from the hash table to FASTQ output files.


## Reference: FASTQ header formats

New R1 header:
- UMI+  : `@A01685:89:HLHWFDRX2:1:1101:4200:1094:GAAAACTC 1:N:0:TTACGGCT+AAGGACCA`
- Note the UMI sequence **`GAAAACTC`**
  
Old R1 header:
- Plain : `@A01685:89:HLHWFDRX2:1:1101:4200:1094 1:N:0:TTACGGCT+AAGGACCA`
- Note: In an R1/R2/R3 system the UMI can technically be present in any read, but it is typically stored in R2.


### Header examples

For UMI **`TCTAAGGC`** and indexes `CTAACTCG+TCGTAGTC`. The actual sequences were clipped only for demonstration purposes.  


#### UMI+ (new system 2021+):

R1  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000:TCTAAGGC 1:N:0:CTAACTCG+TCGTAGTC`
`CNCCAATGTGGAAGTGGATGCTGTAAAATTTAAACTAAAAACACATCTCACCCCAGATGCGTTAGGAGCAAAACGAAGAAGTGCTGGATTCT`  
`+`  
`F#FFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,F:FFFFFFFFFFF,FFFFFFFFFFFFFF:FFFFF,FFFFF`


R2  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000:TCTAAGGC 2:N:0:CTAACTCG+TCGTAGTC`  
`GTAAGAAGTGTCGGTGTATTGGGTGGGTTCGTTCAGATTAAAAATCATTTTAGAATCCAGCACTTCTTCGGTTTGCTCCTAACGCATCTGGG`  
`+`  
`,FFFFF,FFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFF:FF::FFF:F,FF:FFFFF:F`


#### R3 system:

R1  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000 1:N:0:CTAACTCG+TCGTAGTC`
`CNCCAATGTGGAAGTGGATGCTGTAAAATTTAAACTAAAAACACATCTCACCCCAGATGCGTTAGGAGCAAAACGAAGAAGTGCTGGATTCT`  
`+`  
`F#FFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,F:FFFFFFFFFFF,FFFFFFFFFFFFFF:FFFFF,FFFFF`

R2 (UMI)  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000 2:N:0:CTAACTCG+TCGTAGTC`  
**`TCTAAGGC`**  
`+`  
`FFFFFFFF`  

R3  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000 3:N:0:CTAACTCG+TCGTAGTC`  
`GTAAGAAGTGTCGGTGTATTGGGTGGGTTCGTTCAGATTAAAAATCATTTTAGAATCCAGCACTTCTTCGGTTTGCTCCTAACGCATCTGGG`  
`+`  
`,FFFFF,FFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFF:FF::FFF:F,FF:FFFFF:F`  



## Performance


### Processing speed

Processing time is typically limited by disk I/O. The script writes the read pair with the highest summed base quality
from each set of exact duplicates (keeping the best readset).


### Memory usage

Deduplication stores unique reads as keys in memory. Memory usage scales with the number of unique molecules in the dataset.
In practice MAXIMUM memory usage is approximately:
`~2 × size of the uncompressed R1 sequence file`


### Benchmark

Example benchmark on a large dataset:

Dataset: 69 million paired-end reads (clusters) in R1 R2 R3 (R1.gz ≈ 6 GB)

Server: PowerEdge R750, Intel Xeon Gold 6354 (72 threads), 256 GB memory running Ubuntu 24.04.3 LTS.

Results:
- Read + deduplicate: ~6 minutes
- Compression + writing: ~5 minutes
- Total runtime: ~11 minutes

Processing speed is largely limited by gzip (de)compression and disk I/O.
Details in the performance_test directory.


## Author

a.bossers@uu.nl // alex.bossers@wur.nl


## Disclaimer

**Script is provided AS IS under GPL-3.**

We did our best to verify that the results are legitimate. However, the output should be considered erroneous, so you should check your results!
The authors, nor their institutions/employers, are in any way direct or indirect responsible for the direct or indirect damages caused by using this tool.
Use it at your own responsibility.
