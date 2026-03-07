
Seems like the perl script processing these files is as fast as directly reading, decompressing and compressing+writing. So almost no overhead time.
Note on large files: you need to have enough ram to hold the hashtable. his is depeninding on the number of duplicates in the files but roughly calculates to less then 2 x R1 uncompressed.

On a large sample (files R1 R2 R3) containing 96M PE sequences (R1.gz of roughly 7.2GB) it took ~9 minutes to read and dereplicate, and about 14 minutes in addition to gzip and write. Writing helps if you have for instance pigz taking over multithreaded gzip on unix. Which I used.

Total processing time roughly 23 minutes for a sample of 96M PE reads.
See performance folder.

A small DEMO shelscript_runner is also attached for demoing gzip infiles, and allow to manually see that one seq got dereplicated.

Alex 20220915

Update October 2022: Instead of simple overwrite, now the 'best' (highest SUM qualityscore) is kept of the duplicated sequence set and written to disk.
Update August 2023: Script is now compatible with the NEW format of having UMI in the header (R1 and R2 files) instead of older R1 R2 and R3 files.

PERFORMANCE was tested on the initial september version!
