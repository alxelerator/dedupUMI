# Deduplicate FASTQ sequences using UMIs

DedupUMI removes exact duplicate FASTQ read pairs using UMI sequences without requiring alignment to a reference.

Duplicates are identified using concatenation of R1, R2 and UMI sequences.
For duplicate molecules the read pair with the highest total base quality score is retained.

## Why deduplicate?

Exact duplicates (incl UMI) are typically the result of over-amplification / too many amplification cycles in sequence library prepping 'PCR'.  


## Duplicate definition

Reads are considered duplicates when the following combination is identical:

R1 sequence + R2 sequence + UMI

*For the older R3 system this becomes: R1 sequence + R2 sequence + R3 sequence*


## Performance

Processing time is typically limited by disk I/O. The script writes the read pair with the highest summed base quality
from each set of exact duplicates (keeping the best readset).

Latest version allows for input of just R1 and R2 where the UMI should be available in the FASTQ header (new format). 
The older R3 system is still supported as well.  
This R2 or R3 system is handled automatically depending if two or three Rx files are provided respectively.


### Memory usage

Deduplication stores unique read keys in memory.
Memory usage scales with the number of unique molecules in the dataset.

In practice MAXIMUM memory usage is approximately:

- 2 × size of the R1 sequence file
- + optional R3 UMI file


### Benchmark

Example benchmark on a large dataset:

Dataset:

- 96 million paired-end reads (R1.gz ≈ 7.2 GB)

Results:

- Read + deduplicate: ~9 minutes
- Compression + writing: ~14 minutes
- Total runtime: ~23 minutes

Processing speed is therefore largely limited by gzip compression and disk I/O.

Using pigz (parallel gzip) like I did can improve compression speed on multi-core systems.

Additional screenshots can be found in the performance/ folder of this repository.


## FASTQ header formats

New R1 header:
- UMI+  : `@A01685:89:HLHWFDRX2:1:1101:4200:1094:GAAAACTC 1:N:0:TTACGGCT+AAGGACCA`
- Note the UMI sequence **`GAAAACTC`**
  
Old R1 header:
- Plain : `@A01685:89:HLHWFDRX2:1:1101:4200:1094 1:N:0:TTACGGCT+AAGGACCA`
- Note: In an R1/R2/R3 system the UMI can technically be present in any read, but it is typically stored in R2.


## Header examples

For UMI **`TCTAAGGC`** and indexes `CTAACTCG+TCGTAGTC`  


### UMI+ (new system 2021+):

R1  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000:TCTAAGGC 1:N:0:CTAACTCG+TCGTAGTC`
`CNCCAATGTGGAAGTGGATGCTGTAAAATTTAAACTAAAAACACATCTCACCCCAGATGCGTTAGGAGCAAAA.. ..ACTGTAATTGTATCGCCAAAAGCCGAAGAAGTGCTGGATTCT`  
`+`  
`F#FFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,F:FFFFFFFFFFF,FFFFFFF.. ..FFFFFFFFFFFFFF:FFFFFFF:FFFFFFF:FFFFF,FFFFF`


R2  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000:TCTAAGGC 2:N:0:CTAACTCG+TCGTAGTC`  
`GTAAGAAGTGTCGGTGTATTGGGTGGGTTCGTTCAGATTAAAAATCATTTTAGAATCCAGCACTTCTTCGG.. ..ACTGCAGGAGTTTGAATGTATCATTTTGCTCCTAACGCATCTGGG`  
`+`  
`,FFFFF,FFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFF.. ..FFFFFFFFFF,:FFFFF:FFF:,F:FF::FFF:F,FF:FFFFF:F`


### R3 system:

R1  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000 1:N:0:CTAACTCG+TCGTAGTC`
`CNCCAATGTGGAAGTGGATGCTGTAAAATTTAAACTAAAAACACATCTCACCCCAGATGCGTTAGGAGCAAAA.. ..ACTGTAATTGTATCGCCAAAAGCCGAAGAAGTGCTGGATTCT`  
`+`  
`F#FFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,F:FFFFFFFFFFF,FFFFFFF.. ..FFFFFFFFFFFFFF:FFFFFFF:FFFFFFF:FFFFF,FFFFF`

R2 (UMI)  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000 2:N:0:CTAACTCG+TCGTAGTC`  
**`TCTAAGGC`**  
`+`  
`FFFFFFFF`  

R3  

`@A00379:673:HFL2HDRX2:2:1101:1488:1000 3:N:0:CTAACTCG+TCGTAGTC`  
`GTAAGAAGTGTCGGTGTATTGGGTGGGTTCGTTCAGATTAAAAATCATTTTAGAATCCAGCACTTCTTCGG.. ..ACTGCAGGAGTTTGAATGTATCATTTTGCTCCTAACGCATCTGGG`  
`+`  
`,FFFFF,FFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFF.. ..FFFFFFFFFF,:FFFFF:FFF:,F:FF::FFF:F,FF:FFFFF:F`  


## Input and output

NEW system R1 R2 UMI in header:
- Input: R1 and R2 FASTQ readset. UMI should be present in the header of at least R1!
- Output: R1 and R2 filtered for EXACT duplicates based on concatenated seq of R1 R2 keeping highest TOTAL qual score.

Older system R1 R2 R3 UMI in R2:
- Input: R1 R2 and R3 FASTQ readset. For the best total quality to work the UMI should be present in R2!
- Output: R1 R2 R3 filtered for EXACT duplicates based on concatenated seq of R1 R2 R3 keeping highest TOTAL qual score.


## Requirements

- Unix system with `zcat` to allow reading/writing of gzipped files
- FASTQ files in fixed 4-line format having read-pairs in separate files (do not use interleaved files)

Optional but recommended:

- pigz (parallel gzip) for faster compression


## How it works

DedupUMI is implemented as a lightweight Perl script using an in-memory
hash table keyed by the concatenation of:

`R1_sequence + R2_sequence + UMI`

This allows duplicate detection without alignment, sorting, or external
tools.

Because the algorithm only performs sequential FASTQ reading combined
with hash lookups, runtime is typically limited by disk I/O
(decompression, compression and writing FASTQ files), not by the
deduplication logic itself.

For each read pair the script:

1. Extracts the UMI (from the FASTQ header or the R3 file)
2. Constructs a (hash)key: R1_sequence + R2_sequence + UMI
3. If the key is new → the sequence reads and its quality values are stored
4. If the key already exists → the stored sequence reads and its quality scores are replaced only if the
   new readset has a higher total base quality score

After all reads are processed the remaining unique reads are written from the hash table
back to FASTQ output files.


## Important notes

- It writes out the last sequences found of a duplicate set having the highest TOTAL qualityscore.
- Input FASTQ files must follow the standard 4-line FASTQ format, starting at the first record.
- Sequences in R1 R2 (and R3) should be in same order and NOT INTERLEAVED!!


## Author

a.bossers@uu.nl // alex.bossers@wur.nl


## Disclaimer

**Script is provided AS IS under GPL-3.**

We did our best to verify that the results are legitimate. However, the output should be considered erroneous, so you should check your results.
The authors nor their institutions/employers are in any way direct or indirect responsible for the direct or indirect damages caused by using this script.
Use at your own responsibility.
