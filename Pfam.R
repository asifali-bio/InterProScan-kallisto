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

#set e-value for Pfam protein domain match
evalue = 0.05

#select Pfam protein domain annotations less than specified e-value
filteredannotationfiles <- list()
for (i in seq(1:numberofspecies)) {
  filteredannotation = annotationfiles[[i]][which(annotationfiles[[i]]$V9 < evalue),]
  filteredannotationfiles[specieslist[i,]] <- list(filteredannotation)
  rm(filteredannotation)
}

#cycle species by changing the number within double brackets
filteredannotationfiles[[1]]


for (i in seq(1:numberofspecies)) {
  
  #extract Pfam protein domain annotations filtered by e-value
  justGeneP<-filteredannotationfiles[[i]][,c(1,5,6)]
  justGeneTPM<-abundancefiles[[i]][,c(1,5)]
  
  colnames(justGeneP) <- c("gene_id","pfam","domain")
  justGeneP$pfam = as.factor(justGeneP$pfam)
  
  #remove tail end of transcript label
  a<-gsub("(.*)_.*","\\1",justGeneP$gene_id)
  justGeneP$gene_id <- a
  
  colnames(justGeneP)[colnames(justGeneP)=="gene_id"] <- "target_id"
  Data2 = merge(justGeneP, justGeneTPM)
  #just Pfam and TPM
  Data2 <- Data2[c(2,3,4)]
  #sum TPM values
  Data2 = ddply(Data2, c("pfam","domain"), numcolwise(sum))

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
  rm(a, justGeneP, justGeneTPM, Data2, Data3)
}

#table of pooled Pfam protein domains per species filtered by e-value
save(new, new2, file = "Pfam.RData")
load("Pfam.RData")

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
  lapply(names(new2[,-c(1,2)]), \(x)
         filter(new2, if_all(setdiff(names(new2[,-c(1,2)]), x), ~is.na(.)))),
  names(new2[,-c(1,2)]))

for (i in seq(1:numberofspecies)) {
  #remove NA columns
  uniquepfam[[i]] <- uniquepfam[[i]][,c(1,2,2+i)]
}

#unique protein domains per species filtered by e-value with pooled TPM
#cycle species by changing the number within double brackets
uniquepfam[[1]]

#to generate output files
ssannotations <- list()
ssisoforms <- list()
sstranscripts <- list()

for (i in seq(1:numberofspecies)) {
  
  trinity = c()
  annotation = c()
  
  for (j in seq(1:dim(uniquepfam[[i]])[1])) {
    a1 = uniquepfam[[i]][j,1]
    #grab first unique protein domain entry
    a2 = grepl(a1, filteredannotationfiles[[i]]$V5)
    #locate position of unique protein domain in filtered data
    a3 = which(a2)
    #extract transcript ID
    a4 = as.character(filteredannotationfiles[[i]][a3,1])
    #extract Pfam protein domain ID
    a5 = as.character(filteredannotationfiles[[i]][a3,5])
    #build data frame of transcripts with unique protein domains filtered by e-value
    trinity = c(trinity, a4)
    #build data frame of unique protein domains filtered by e-value
    annotation = c(annotation, a5)
  }
  
  annotation = as.data.frame(annotation)
  #label species
  colnames(annotation) = specieslist[i,]
  #unique protein domains
  annotation = unique(annotation)

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
  
  
  rm(a1, a2, a3, a4, a5, clean, clean2, annotation, trinity, trinity2)
}


#species-specific information
#cycle species by changing the number within double brackets

#unique protein domains filtered by e-value
ssannotations[[1]]

#transcript isoforms with unique protein domains filtered by e-value
ssisoforms[[1]]

#source genes of transcript isoforms with unique protein domains filtered by e-value
sstranscripts[[1]]


#set working directory to new folder

#save
for (i in seq(1:numberofspecies)) {
  #species-specific Pfam protein domains
  write.table(uniquepfam[[i]][,c(1,2)], file = paste0(specieslist[i,], "_Pfam.txt"), col.names = FALSE)
  #species-specific transcript isoforms
  write.table(ssisoforms[[i]], file = paste0(specieslist[i,], "_i.txt"), col.names = FALSE)
  #species-specific genes
  write.table(sstranscripts[[i]], file = paste0(specieslist[i,], "_g.txt"), col.names = FALSE)
}



#Part C
#evolutionary distance based on clustering

#1
#cluster by species
new3 = new2
rownames(new3) = new3[,1]
new3 = new3[,-c(1,2)]
pheatmap(new3, scale = "row", treeheight_row = 0, cluster_cols = TRUE, cluster_rows = FALSE, show_rownames = FALSE, na_col = "black")

#2
#cluster by species and annotation
new3 = new2
rownames(new3) = new3[,1]
new3 = new3[,-c(1,2)]
#only complete cases
new3 = new3[complete.cases(new3), ]
pheatmap(new3, scale = "row", treeheight_row = 0, cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = FALSE)

#3
#cluster by species and annotation
new3 = new2
rownames(new3) = new3[,1]
new3 = new3[,-c(1,2)]
#set NA to 0
new3[is.na(new3)] = 0
#calculate distance matrix
d = dist(new3, method = "euclidean")
#cluster
hc = hclust(d, method = "complete")
order = hc$labels[hc$order]
#rearrange data
new4 = new2 %>%
  arrange(factor(pfam, levels = order))
rownames(new4) = new4[,1]
new4 = new4[,-c(1,2)]
#plot NA values as NA
pheatmap(new4, scale = "row", treeheight_row = 0, cluster_cols = TRUE, cluster_rows = FALSE, show_rownames = FALSE, na_col = "black")
