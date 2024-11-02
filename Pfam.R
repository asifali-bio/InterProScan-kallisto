library(plyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(viridis)


#read CSV file containing list of sample/species names
specieslist = read.csv("specieslist.csv", header = F)
numberofspecies = nrow(specieslist)

#read kallisto output files as a list
abundancelist = paste0("abundance", 1:numberofspecies, ".tsv")
abundancefiles = lapply(abundancelist, read.delim)

#read InterProScan output files as a list
annotationlist = paste0("annotation", 1:numberofspecies, ".tsv")
annotationfiles = lapply(annotationlist, read.delim, header=F)

#set e-value for Pfam protein domain match
evalue = 0

#remove Pfam protein domain annotations less than or equal to e-value
filteredannotationfiles <- list()
for (i in seq(1:numberofspecies)) {
  filteredannotation = annotationfiles[[i]][which(annotationfiles[[i]]$V9 >= evalue),]
  filteredannotationfiles[specieslist[i,]] <- list(filteredannotation)
}

filteredannotationfiles[[1]]
#cycle species by changing the number within double brackets


for (i in seq(1:numberofspecies)) {
  
  #extract Pfam protein domain annotations filtered by e-value
  justGeneP<-filteredannotationfiles[[i]][,c(1,5)]
  justGeneTPM<-abundancefiles[[i]][,c(1,5)]
  
  colnames(justGeneP) <- c("gene_id","pfam")
  justGeneP$pfam = as.factor(justGeneP$pfam)
  
  #remove tail end of transcript label
  a<-gsub("(.*)_.*","\\1",justGeneP$gene_id)
  justGeneP$gene_id <- a
  
  colnames(justGeneP)[colnames(justGeneP)=="gene_id"] <- "target_id"
  Data2 = merge(justGeneP, justGeneTPM)
  Data2 <- Data2[c(2,3)]
  #just Pfam and TPM
  Data2 = ddply(Data2, "pfam", numcolwise(sum))
  #sum TPM values
  
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
  rm(a, justGeneP, justGeneTPM, Data2, Data3)
}

#table of pooled Pfam protein domains per species filtered by e-value
save(new, new2, file = "Pfam.RData")
load("Pfam.RData")

#calculate standard score for heatmap
#new3 = t(scale(t(new2[,-1])))
#new3[is.nan(new3)] <- 0
#new3 = as.list(new3)
#new3[is.na(new3)] <- NULL
#for (i in seq(1:nrow(new))) {
#  new[i,2] = new3[i]
#}

#preview plot
ggplot(new, aes(species, pfam)) +
  geom_point(aes(color = pfam, size = tpm), alpha = 0.3, show.legend = FALSE) +
  theme_classic() +
  labs(title = "Derived Pfam protein domain distributions across species", x = "Species", y = "Pfam") +
  scale_fill_gradientn(colours = rainbow(nrow(new2))) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

#plot of pooled Pfam protein domains per species filtered by e-value
p1 = plot_ly(x = new$species, y = new$pfam, type = "scatter", mode = "markers", color = new$species, size = new$tpm, fill = ~'', colors = "Spectral")
p1 <- p1 %>% layout(title = "InterProScan x kallisto",
                    xaxis = list(title = "Species"),
                    yaxis = list(title = "Pfam", showticklabels = FALSE),
                    showlegend = FALSE)
#color by species
p1


p2 = plot_ly(x = new$species, y = new$pfam, type = "scatter", mode = "markers", color = new$species, size = new$tpm, fill = ~'', colors = "Spectral")
p2 <- p2 %>% layout(title = "InterProScan x kallisto",
                    xaxis = list(title = "Species", showticklabels = FALSE),
                    yaxis = list(title = "Pfam", showticklabels = FALSE))
#color by species
p2


p3 = plot_ly(x = new$species, y = new$pfam, type = "scatter", mode = "markers", color = new$pfam, size = new$tpm, fill = ~'', colors = viridis(nrow(new2), direction = -1))
p3 <- p3 %>% layout(title = "InterProScan x kallisto",
                    xaxis = list(title = "Species"),
                    yaxis = list(title = "Pfam", showticklabels = FALSE),
                    showlegend = FALSE)
#color by annotation
p3



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


ssannotations <- list()
ssisoforms <- list()
sstranscripts <- list()

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
    #locate position of unique protein domain in filtered data
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
  annotation = unique(annotation)
  #unique protein domains
  colnames(annotation) = specieslist[i,]
  #label species
  
  trinity = as.data.frame(trinity)
  trinity2 = trinity
  
  clean <- gsub("(.*)_.*","\\1",trinity$trinity)
  trinity$trinity <- clean
  #remove tail end of transcript label
  trinity = unique(trinity)
  #unique transcripts
  colnames(trinity) = specieslist[i,]
  #label species
  
  clean <- gsub("(.*)_.*","\\1",trinity2$trinity)
  trinity2$trinity <- clean
  clean2 <- gsub("(.*)_.*","\\1",trinity2$trinity)
  #remove isoform tag
  trinity2$trinity <- clean2
  trinity2 = unique(trinity2)
  colnames(trinity2) = specieslist[i,]
  
  #create a list of annotations and transcripts for each species
  ssannotations[specieslist[i,]] <- list(annotation)
  ssisoforms[specieslist[i,]] <- list(trinity)
  sstranscripts[specieslist[i,]] <- list(trinity2)
  
  
  rm(a, a1, a2, a3, a4, a5, clean, clean2, filteredannotation, annotation, trinity, trinity2)
}


#species-specific information
#cycle species by changing the number within double brackets

ssannotations[[1]]
#unique protein domains filtered by e-value

ssisoforms[[1]]
#transcript isoforms with unique protein domains filtered by e-value

sstranscripts[[1]]
#source genes of transcript isoforms with unique protein domains filtered by e-value

#set working directory to new folder
for (i in seq(1:numberofspecies)) {
  write.table(ssannotations[[i]], file = paste0(specieslist[i,], "_Pfam.txt"), col.names = FALSE)
  write.table(ssisoforms[[i]], file = paste0(specieslist[i,], "_i.txt"), col.names = FALSE)
  write.table(sstranscripts[[i]], file = paste0(specieslist[i,], "_g.txt"), col.names = FALSE)
  #save
}
#species-specific Pfam protein domains
#species-specific transcript isoforms
#species-specific genes