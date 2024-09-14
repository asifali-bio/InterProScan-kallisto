#install.packages("BiocManager")
#BiocManager::install("GOstats")

#library("GOstats")
#library("GSEABase")
#library(GO.db)
library(ggplot2)
library(plotly)
library(plyr)

annotation = read.delim("ABIJ.tsv", header=F)
head(annotation)
annotation2 = read.delim("JKAA.tsv", header=F)
annotation3 = read.delim("KJYC.tsv", header=F)
annotation4 = read.delim("KUXM.tsv", header=F)
annotation5 = read.delim("LGDQ.tsv", header=F)
annotation6 = read.delim("ZFGK.tsv", header=F)
annotation7 = read.delim("ZYCD.tsv", header=F)
annotation8 = read.delim("ZZOL.tsv", header=F)


justGeneGo<-annotation[,c(1,14)]

go_id=character()
gene_id=character()
#empty character vectors

#loop to collect all go ids associated with each scaffold
simple_counter=1
for (i in 1:nrow(justGeneGo)) {
  gene=as.character(justGeneGo[i,1])
  
  
  #multiple GO-terms are delimited "|"
  go=unlist(strsplit(as.character(justGeneGo[i,2]),"|", fixed=T))
  
  
  for (eachgo in go) {
    for (eachgene in gene) {
      go_id[simple_counter]=eachgo
      gene_id[simple_counter]=eachgene
      simple_counter=simple_counter+1
    }
  }
}

Data = data.frame(go_id, gene_id)


a<-gsub("(.*)_.*","\\1",Data$gene_id)
Data$gene_id <- a


#justGeneP<-annotation[,c(1,5)]
justGeneTPM<-ABIJ.abundance[,c(1,5)]


colnames(Data)[colnames(Data)=="gene_id"] <- "target_id"
Data2 = merge(Data, justGeneTPM)
#Data2 <- Data2[c(2,3)]
#just GO and TPM
#Data2 = ddply(Data2, "go_id", numcolwise(sum))
#sum tpm values

#Data3 = Data2

#colnames(Data3)[colnames(Data3)=="tpm"] <- "lepidophylla"
#Data2$newcolumn <- "lepidophylla"
#label


#Selaginella lepidophylla (ABIJ)
#Selaginella wallacei (JKAA)
#Selaginella willdenowii (KJYC)
#Selaginella selaginoides (KUXM)
#Selaginella apoda (LGDQ)
#Selaginella kraussiana (ZFGK)
#Selaginella acanthonota (ZYCD)
#Selaginella stauntoniana (ZZOL)



#Data$go_id = as.character(Data$go_id)
#Data$go_id = Term(Data$go_id)
#Data$go_id = as.factor(Data$go_id)

#plot(Data$go_id, ylim=c(0, 12000))


#goterms <- Term(GOTERM)
#head(goterms)

#goterms["GO:0006468"]

#GOID("GO:0006468")
#Term("GO:0006468")
#Synonym("GO:0006468")
#Secondary("GO:0006468")
#Definition("GO:0006468")
#Ontology("GO:0006468")


new2 = Data3
#new2 = merge.data.frame(new2, Data3, all = TRUE)


new = Data2
#new = rbind(new, Data2)
#combine all data

#p
new$go_id = as.character(new$go_id)
#new$go_id = Ontology(new$go_id)

#new = Data
#add new column for count per GO term
#new3 = ddply(new,.(go_id, gene_id), nrow)
#colnames(new3)[3]<-"count"


#new3$go_id = as.factor(new3$go_id)
#new3$gene_id = as.factor(new3$gene_id)


p = plot_ly(x = new$newcolumn, y = new$go_id, type = "scatter", mode = "markers", color = new$go_id, size = new$tpm, colors = "Spectral")
p


p1 = plot_ly(x = new$go_id, y = new$tpm, type = "scatter", mode = "markers", color = new$tpm, size = new$tpm, colors = "Blues")
p1


save(new, new2, file = "all.RData")

load("all.RData")


a = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
b = new2[is.na(new2$lepidophylla) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
c = new2[is.na(new2$wallacei) & is.na(new2$lepidophylla) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
d = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$lepidophylla) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
e = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$lepidophylla) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
f = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$lepidophylla) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
g = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$lepidophylla) & is.na(new2$stauntoniana),]
h = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$lepidophylla),]



new2$go_id <- NULL
colSums(!is.na(new2))

new2[is.na(new2)] <- 0
colSums(new2)

new2 = new2[complete.cases(new2), ]



trinity = c()
annot = c()

a5 = dim(h)[1]

for (i in seq(1:a5)) {
  a3 = h[i,1]
  a1 = grepl(a3, annotation8$V14)
  a2 = which(a1)
  a4 = as.character(annotation8[a2,1])
  a5 = as.character(annotation8[a2,14])
  trinity = c(trinity, a4)
  annot = c(annot, a5)
}

annot = as.data.frame(annot)

trinity = as.data.frame(trinity)
clean <- gsub("(.*)_.*","\\1",trinity$trinity)
trinity$trinity <- clean
#clean2 <- gsub("(.*)_.*","\\1",trinity$trinity)
#trinity$trinity <- clean2
trinity = unique(trinity)

write.table(trinity, file = "a8.txt", col.names = FALSE)
