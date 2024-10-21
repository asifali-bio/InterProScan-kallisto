library(ggplot2)
library(plotly)
library(plyr)
library(dplyr)


#make a list of sample names in CSV format via Excel
specieslist = read.csv("specieslist.csv", header = F)
numberofspecies = nrow(specieslist)

abundancelist = paste0("abundance", 1:numberofspecies, ".tsv")
abundancefiles = lapply(abundancelist, read.delim)
#abundance1 = as.data.frame(abundancefiles[1])

annotationlist = paste0("annotation", 1:numberofspecies, ".tsv")
annotationfiles = lapply(annotationlist, read.delim, header=F)
#annotation1 = as.data.frame(annotationfiles[1])


for (i in seq(1:numberofspecies)) {
  
  justGeneP<-as.data.frame(annotationfiles[i])[,c(1,5)]
  justGeneTPM<-as.data.frame(abundancefiles[i])[,c(1,5)]
  
  colnames(justGeneP) <- c("gene_id","pfam") 
  
  a<-gsub("(.*)_.*","\\1",justGeneP$gene_id)
  justGeneP$gene_id <- a
  
  colnames(justGeneP)[colnames(justGeneP)=="gene_id"] <- "target_id"
  Data2 = merge(justGeneP, justGeneTPM)
  Data2 <- Data2[c(2,3)]
  #just pfam and TPM
  Data2 = ddply(Data2, "pfam", numcolwise(sum))
  #sum tpm values
  
  Data3 = Data2
  
  colnames(Data3)[colnames(Data3)=="tpm"] <- specieslist[i,]
  Data2$newcolumn <- specieslist[i,]
  #label
  
  if (i==1) {
    new2 = Data3
    new = Data2
    #initialize data
  }
  else {
    new2 = merge.data.frame(new2, Data3, all = TRUE)
    new = rbind(new, Data2)
    #combine all data
  }
}


p1 = plot_ly(x = new$newcolumn, y = new$pfam, type = "scatter", mode = "markers", color = new$newcolumn, size = new$tpm, fill = ~'', colors = "Spectral")
p1 <- p1 %>% layout(showlegend = FALSE)

p1



#Part B

uniquepfam = setNames(
  lapply(names(new2[,-1]), \(x)
         filter(new2, if_all(setdiff(names(new2[,-1]), x), ~is.na(.)))),
  names(new2[,-1]))

#a=as.data.frame(uniquepfam[1])
#unique protein domains

#annotation = as.data.frame(annotationfiles[1])[which(as.data.frame(annotationfiles[1])$V9 < 0.001),]
#filter by e-value

finalannotationlist <- list()
finaltranscriptlist <- list()

for (i in seq(1:numberofspecies)) {
  a=as.data.frame(uniquepfam[i])
  annotation = as.data.frame(annotationfiles[i])[which(as.data.frame(annotationfiles[i])$V9 < 0.001),]
  
  trinity = c()
  annot = c()
  
  a5 = dim(a)[1]
  #change letter according to sample
  
  for (j in seq(1:a5)) {
    a3 = a[j,1]
    #grab first unique protein domain entry
    a1 = grepl(a3, annotation$V5)
    a2 = which(a1)
    #grab position of unique protein domain in filtered data
    a4 = as.character(annotation[a2,1])
    #extract transcript ID
    a5 = as.character(annotation[a2,5])
    #extract Pfam protein domain ID
    trinity = c(trinity, a4)
    #build data frame of transcripts with unique protein domains filtered by E-value
    annot = c(annot, a5)
    #build data frame of unique protein domains filtered by E-value
  }
  
  annot = as.data.frame(annot)
  annot2 = unique(annot)
  #unique protein domains
  colnames(annot2) = specieslist[i,]
  #label species
  
  trinity = as.data.frame(trinity)
  clean <- gsub("(.*)_.*","\\1",trinity$trinity)
  trinity$trinity <- clean
  #remove tail end of transcript label
  trinity2 = unique(trinity)
  #unique transcripts
  colnames(trinity2) = specieslist[i,]
  #label species
  
  #create a list of annotations and transcripts for each species
  finalannotationlist[specieslist[i,]] <- list(annot2)
  finaltranscriptlist[specieslist[i,]] <- list(trinity2)
}


#a = as.data.frame(finalannotationlist[1])
#b = as.data.frame(finaltranscriptlist[1])



