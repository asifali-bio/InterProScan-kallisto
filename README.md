For de novo RNA-seq analysis, Trinity was used to assemble reads. Then, we applied InterProScan as well as kallisto to each Trinity file. Once we have the output files from InterProScan and kallisto, we combine the GO/Pfam factors with abundance information in R. The R code contained in this repository assumes that both InterProScan and kallisto output files are ready to be processed and used as input files.

In theory, we should be able to perform RNA-seq analysis across different species. This method treats evolutionary time as a treatment condition. As such, there is no normalization step. Normalization would assume that most homologs are not differentially expressed and that total gene expression is equal across species. Keep in mind that read counts measure relative abundance per sample, not absolute abundance.

The main purpose of this method is to identify transcripts with a unique protein domain signature across different species.

An analogue sequence similarity algorithm has been developed and will be implemented in future updates.


To run the function, name the abundance and annotation input files in a similar format as well as the species list CSV file. In R, set the working directory to the same location as the input files. Run the code line-by-line.