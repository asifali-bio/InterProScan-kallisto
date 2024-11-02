library(plyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(viridis)


#read CSV containing list of sample/species names
specieslist = read.csv("specieslist.csv", header = F)
numberofspecies = nrow(specieslist)

#read kallisto output files as a list
abundancelist = paste0("abundance", 1:numberofspecies, ".tsv")
abundancefiles = lapply(abundancelist, read.delim)

#read InterProScan output files as a list
annotationlist = paste0("annotation", 1:numberofspecies, ".tsv")
annotationfiles = lapply(annotationlist, read.delim, header=F)


for (i in seq(1:numberofspecies)) {
  
  justGeneGo<-annotationfiles[[i]][,c(1,14)]
  justGeneTPM<-abundancefiles[[i]][,c(1,5)]
  
  #empty character vectors
  go_id=character()
  gene_id=character()

  #loop to collect all GO terms associated with each scaffold
  simple_counter=1
  for (j in 1:nrow(justGeneGo)) {
    gene=as.character(justGeneGo[j,1])
    
    #multiple GO terms are delimited "|"
    go=unlist(strsplit(as.character(justGeneGo[j,2]),"|", fixed=T))
    
    
    for (eachgo in go) {
      for (eachgene in gene) {
        go_id[simple_counter]=eachgo
        gene_id[simple_counter]=eachgene
        simple_counter=simple_counter+1
      }
    }
  }
  
  go_id = as.factor(go_id)
  Data = data.frame(go_id, gene_id)
  
  #remove tail end of transcript label
  a<-gsub("(.*)_.*","\\1",Data$gene_id)
  Data$gene_id <- a
  
  colnames(Data)[colnames(Data)=="gene_id"] <- "target_id"
  Data2 = merge(Data, justGeneTPM)
  Data2 <- Data2[c(2,3)]
  #just GO and TPM
  Data2 = ddply(Data2, "go_id", numcolwise(sum))
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
  rm(a, eachgene, eachgo, gene, go, gene_id, go_id, simple_counter, justGeneGo, justGeneTPM, Data, Data2, Data3)
}

#table of pooled GO terms per species
save(new, new2, file = "GO.RData")

load("GO.RData")
ggplot(new, aes(species, go_id)) +
  geom_point(aes(color = go_id, size = tpm), alpha = 0.3, show.legend = FALSE) +
  theme_classic() +
  labs(title = "Derived GO term distributions across species", x = "Species", y = "GO") +
  scale_fill_gradientn(colours = rainbow(nrow(new2))) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

#plot of pooled GO terms per species
p1 = plot_ly(x = new$species, y = new$go_id, type = "scatter", mode = "markers", color = new$species, size = new$tpm, fill = ~'', colors = "Spectral")
p1 <- p1 %>% layout(title = "InterProScan x kallisto",
                    xaxis = list(title = "Species"),
                    yaxis = list(title = "GO", showticklabels = FALSE),
                    showlegend = FALSE)
#color by species
p1


p2 = plot_ly(x = new$species, y = new$go_id, type = "scatter", mode = "markers", color = new$species, size = new$tpm, fill = ~'', colors = "Spectral")
p2 <- p2 %>% layout(title = "InterProScan x kallisto",
                    xaxis = list(title = "Species", showticklabels = FALSE),
                    yaxis = list(title = "GO", showticklabels = FALSE))
#color by species
p2


p3 = plot_ly(x = new$species, y = new$go_id, type = "scatter", mode = "markers", color = new$go_id, size = new$tpm, fill = ~'', colors = viridis(nrow(new2), direction = -1))
p3 <- p3 %>% layout(title = "InterProScan x kallisto",
                    xaxis = list(title = "Species"),
                    yaxis = list(title = "GO", showticklabels = FALSE),
                    showlegend = FALSE)
#color by annotation
p3



#Part B

uniqueGO = setNames(
  lapply(names(new2[,-1]), \(x)
         filter(new2, if_all(setdiff(names(new2[,-1]), x), ~is.na(.)))),
  names(new2[,-1]))

#unique GO terms per species
#corresponding pooled TPM included
uniqueGO[[1]]
#cycle species by changing the number within double brackets


ssannotations <- list()
ssisoforms <- list()
sstranscripts <- list()

for (i in seq(1:numberofspecies)) {
  a=uniqueGO[[i]]
  annotationGo = annotationfiles[[i]]
  
  trinity = c()
  annotation = c()
  
  a5 = dim(a)[1]
  #set length of loop
  
  for (j in seq(1:a5)) {
    a3 = a[j,1]
    #grab first unique GO term entry
    a1 = grepl(a3, annotationGo$V14)
    a2 = which(a1)
    #locate position of unique GO term
    a4 = as.character(annotationGo[a2,1])
    #extract transcript ID
    trinity = c(trinity, a4)
    #build data frame of transcripts with unique GO terms
  }
  
  annotation = a[,1]
  annotation = as.data.frame(annotation)
  annotation = unique(annotation)
  #unique GO terms
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
  
  
  rm(a, a1, a2, a3, a4, a5, clean, clean2, annotationGo, annotation, trinity, trinity2)
}


#species-specific information
#cycle species by changing the number within double brackets

ssannotations[[1]]
#unique GO terms

ssisoforms[[1]]
#transcript isoforms with unique GO terms

sstranscripts[[1]]
#source genes of transcript isoforms with unique GO terms

#set working directory to new folder
for (i in seq(1:numberofspecies)) {
  write.table(ssannotations[[i]], file = paste0(specieslist[i,], "_GO.txt"), col.names = FALSE)
  write.table(ssisoforms[[i]], file = paste0(specieslist[i,], "_i.txt"), col.names = FALSE)
  write.table(sstranscripts[[i]], file = paste0(specieslist[i,], "_g.txt"), col.names = FALSE)
  #save
}
#species-specific GO terms
#species-specific transcript isoforms
#species-specific genes