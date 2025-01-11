library(plyr)
library(dplyr)
library(pheatmap)
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



#Part A

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
  #just GO and TPM
  Data2 <- Data2[c(2,3)]
  #sum TPM values
  Data2 = ddply(Data2, "go_id", numcolwise(sum))

  Data3 = Data2
  
  colnames(Data3)[colnames(Data3)=="tpm"] <- specieslist[i,]
  #label
  Data2$species <- specieslist[i,]

  if (i==1) {
    #initialize data
    new2 = Data3
    new = Data2
  }
  else {
    #combine all data
    new2 = merge.data.frame(new2, Data3, all = TRUE)
    new = rbind(new, Data2)
  }
  rm(a, eachgene, eachgo, gene, go, gene_id, go_id, simple_counter, justGeneGo, justGeneTPM, Data, Data2, Data3)
}

#table of pooled GO terms per species
save(new, new2, file = "GO.RData")
load("GO.RData")

#preview plot
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

for (i in seq(1:numberofspecies)) {
  #remove NA columns
  uniqueGO[[i]] <- uniqueGO[[i]][,c(1,1+i)]
}

#unique GO terms per species with pooled TPM
#cycle species by changing the number within double brackets
uniqueGO[[1]]


ssannotations <- list()
ssisoforms <- list()
sstranscripts <- list()

for (i in seq(1:numberofspecies)) {
  
  trinity = c()
  
  for (j in seq(1:dim(uniqueGO[[i]])[1])) {
    a1 = uniqueGO[[i]][j,1]
    #grab first unique GO term entry
    a2 = grepl(a1, annotationfiles[[i]]$V14)
    #locate position of unique GO term
    a3 = which(a2)
    #extract transcript ID
    a4 = as.character(annotationfiles[[i]][a3,1])
    #build data frame of transcripts with unique GO terms
    trinity = c(trinity, a4)
  }
  
  annotation = uniqueGO[[i]][,1]
  annotation = as.data.frame(annotation)
  #label species
  colnames(annotation) = specieslist[i,]

  trinity = as.data.frame(trinity)
  #label species
  colnames(trinity) = specieslist[i,]
  #remove tail end of transcript label
  clean <- gsub("(.*)_.*","\\1",trinity[,1])
  trinity[,1] <- clean
  #unique transcripts
  trinity = unique(trinity)
  
  #make a copy
  trinity2 = trinity
  #remove isoform tag
  clean2 <- gsub("(.*)_.*","\\1",trinity2[,1])
  trinity2[,1] <- clean2
  #unique isoforms
  trinity2 = unique(trinity2)
  
  #create a list of annotations and transcripts for each species
  ssannotations[specieslist[i,]] <- list(annotation)
  ssisoforms[specieslist[i,]] <- list(trinity)
  sstranscripts[specieslist[i,]] <- list(trinity2)
  
  
  rm(a1, a2, a3, a4, clean, clean2, annotation, trinity, trinity2)
}


#species-specific information
#cycle species by changing the number within double brackets

#unique GO terms
ssannotations[[1]]

#transcript isoforms with unique GO terms
ssisoforms[[1]]

#source genes of transcript isoforms with unique GO terms
sstranscripts[[1]]


#set working directory to new folder

#save
for (i in seq(1:numberofspecies)) {
  #species-specific GO terms
  write.table(ssannotations[[i]], file = paste0(specieslist[i,], "_GO.txt"), col.names = FALSE)
  #species-specific transcript isoforms
  write.table(ssisoforms[[i]], file = paste0(specieslist[i,], "_i.txt"), col.names = FALSE)
  #species-specific genes
  write.table(sstranscripts[[i]], file = paste0(specieslist[i,], "_g.txt"), col.names = FALSE)
}



#Part C
#evolutionary distance based on clustering

remove = c()
#unique GO terms are irrelevant for clustering

for (i in seq(1:numberofspecies)) {
  for (j in seq(1:dim(uniqueGO[[i]])[1])) {
    a1 = as.character(uniqueGO[[i]][j,1])
    #build list of unique GO terms for removal
    remove = c(remove, a1)
  }
  rm(a1)
}

#remove unique GO terms
new3 = new2
new3 = new3[!new3$go_id %in% remove,]
new3 = new3[,-1]
pheatmap(new3, scale = "row", treeheight_row = 0, cluster_cols = TRUE, cluster_rows = FALSE, show_rownames = FALSE, na_col = "black")

#only complete cases for clustering species and GO terms
new3 = new2
new3 = new3[complete.cases(new3), ]
new3 = new3[,-1]
pheatmap(new3, scale = "row", treeheight_row = 0, cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = FALSE)

#set NA to 0 for clustering species and GO terms
new3 = new2
new3[is.na(new3)] = 0
new3 = new3[,-1]
pheatmap(new3, scale = "row", treeheight_row = 0, cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = FALSE)
