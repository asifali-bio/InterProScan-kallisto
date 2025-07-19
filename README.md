\*\*Functional genomics and comparative systems biology\*\*



This pipeline allows for de novo RNA-seq analysis across different species using structural annotation and expression quantification. Reads were assembled via Trinity. Then, protein domains were predicted and gene expression was estimated via InterProScan and kallisto, respectively. Once we have the output files from InterProScan and kallisto, we use R to merge the GO/Pfam annotations with the abundance information to obtain a distribution of annotations. The R code contained in this repository assumes that both InterProScan and kallisto output files are ready to be processed and used as input files.



This method treats evolutionary time as a treatment condition. There is no normalization step due to the vast differences in expression profiles across different species. Normalization would assume that most homologs are not differentially expressed and that total gene expression is equal across species. Quantification of the mRNA expression levels allows us to assess within species variance, as read counts measure relative abundance per species rather than absolute abundance.



The main purpose of this method is to identify species-specific transcripts that translate to proteins containing a protein domain not found in other species.



To run the function, name the abundance and annotation input files in a similar format as well as the species list CSV file. In R, set the working directory to the same location as the input files. Run the code line by line.



An analogue sequence similarity algorithm has been developed and will be implemented in a future update.

