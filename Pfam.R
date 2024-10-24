library(ggplot2)
library(plotly)
library(plyr)
library(dplyr)


#read CSV containing list of sample/species names
specieslist = read.csv("specieslist.csv", header = F)
numberofspecies = nrow(specieslist)

#read kallisto output files as a list
abundancelist = paste0("abundance", 1:numberofspecies, ".tsv")
abundancefiles = lapply(abundancelist, read.delim)

#read InterProScan output files as a list
annotationlist = paste0("annotation", 1:numberofspecies, ".tsv")
annotationfiles = lapply(annotationlist, read.delim, header=F)

#set e-value threshold for Pfam protein domain match
evalue = 0.001

#filter annotation files by e-value
filteredannotationfiles <- list()
for (i in seq(1:numberofspecies)) {
  filteredannotation = annotationfiles[[i]][which(annotationfiles[[i]]$V9 < evalue),]
  filteredannotationfiles[specieslist[i,]] <- list(filteredannotation)
}

filteredannotationfiles[[1]]
#cycle species by changing the number within double brackets


for (i in seq(1:numberofspecies)) {
  
  #filter Pfam protein domain annotations by e-value
  justGeneP<-filteredannotationfiles[[i]][,c(1,5)]
  justGeneTPM<-abundancefiles[[i]][,c(1,5)]
  
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
  Data2$species <- specieslist[i,]
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

#table of pooled protein domains per species filtered by e-value
new2

#plot of pooled protein domains per species filtered by e-value
p1 = plot_ly(x = new$species, y = new$pfam, type = "scatter", mode = "markers", color = new$species, size = new$tpm, fill = ~'', colors = "Spectral")
p1 <- p1 %>% layout(showlegend = FALSE)

p1



#Part B

uniquepfam = setNames(
  lapply(names(new2[,-1]), \(x)
         filter(new2, if_all(setdiff(names(new2[,-1]), x), ~is.na(.)))),
  names(new2[,-1]))

#unique protein domains per species
#corresponding pooled TPM included
#filtered by e-value
uniquepfam[[1]]
#cycle species by changing the number within double brackets


finalannotationlist <- list()
finaltranscriptlist <- list()

for (i in seq(1:numberofspecies)) {
  a=uniquepfam[[i]]
  filteredannotation = filteredannotationfiles[[i]]
  
  trinity = c()
  annotation = c()
  
  a5 = dim(a)[1]
  #set length of loop
  
  for (j in seq(1:a5)) {
    a3 = a[j,1]
    #grab first unique protein domain entry
    a1 = grepl(a3, filteredannotation$V5)
    a2 = which(a1)
    #grab position of unique protein domain in filtered data
    a4 = as.character(filteredannotation[a2,1])
    #extract transcript ID
    a5 = as.character(filteredannotation[a2,5])
    #extract Pfam protein domain ID
    trinity = c(trinity, a4)
    #build data frame of transcripts with unique protein domains filtered by e-value
    annotation = c(annotation, a5)
    #build data frame of unique protein domains filtered by e-value
  }
  
  annotation = as.data.frame(annotation)
  annotation2 = unique(annotation)
  #unique protein domains
  colnames(annotation2) = specieslist[i,]
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
  finalannotationlist[specieslist[i,]] <- list(annotation2)
  finaltranscriptlist[specieslist[i,]] <- list(trinity2)
}


#access species-specific information
#cycle species by changing the number within double brackets

finalannotationlist[[1]]
#unique protein domains filtered by e-value

finaltranscriptlist[[1]]
#transcripts with unique protein domains filtered by e-value

