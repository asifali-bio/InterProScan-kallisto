1: On Unix, download raw RNA-seq reads from SRA database, each sample split into two independent FASTQ files. Run FastQC on raw RNA-seq reads.

2: Assemble transcriptome via Trinity on HPC.

3: Calculate gene expression (read counts) via kallisto on HPC.

4a: On Unix, prepare data for InterProScan parallel processing. Trim header. We want to split the single FASTA file into many FASTA files that can be processed in parallel, so we create a folder where each FASTA file only contains 100 sequences. To lower the thread count for decreased parallelism, we can increase the number of sequences per FASTA file.

4b: Annotate transcriptome via InterProScan parallel processing on HPC.

4c: On Unix, merge all TSV files created by parallel processing job after navigating to a single folder.