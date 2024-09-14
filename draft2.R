#install.packages("BiocManager")
#BiocManager::install("GOstats")

library("GOstats")
library("GSEABase")
library(GO.db)
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

annotation <- annotation[which(annotation$V9 < 0.001),]
annotation2 <- annotation2[which(annotation2$V9 < 0.001),]
annotation3 <- annotation3[which(annotation3$V9 < 0.001),]
annotation4 <- annotation4[which(annotation4$V9 < 0.001),]
annotation5 <- annotation5[which(annotation5$V9 < 0.001),]
annotation6 <- annotation6[which(annotation6$V9 < 0.001),]
annotation7 <- annotation7[which(annotation7$V9 < 0.001),]
annotation8 <- annotation8[which(annotation8$V9 < 0.001),]


justGeneP<-annotation8[,c(1,5,9)]
#load kallisto abundance TSV
justGeneTPM<-ZZOL.abundance[,c(1,5)]

colnames(justGeneP) <- c("gene_id","pfam","e_value") 

a<-gsub("(.*)_.*","\\1",justGeneP$gene_id)
justGeneP$gene_id <- a
a = 0

colnames(justGeneP)[colnames(justGeneP)=="gene_id"] <- "target_id"
Data2 = merge(justGeneP, justGeneTPM)
#Data2 <- Data2[c(2,3,4)]

#just pfam and TPM
#Data2 = ddply(Data2, "pfam", numcolwise(sum))
#sum tpm values

#Data3 = Data2

#colnames(Data3)[colnames(Data3)=="tpm"] <- "stauntoniana"
#Data2$newcolumn <- "stauntoniana"
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


#new2 = Data3
#new2 = merge.data.frame(new2, Data3, all = TRUE)


new = Data2
#new = rbind(new, Data2)
#combine all data


#Data2$count <- 1
#Data2 <- Data2[order(Data2$e_value),]

#track2 = c()
#prot = c()
#for (i in Data2$pfam) {
#  if (is.element(i, prot) == TRUE) {
#    track = sum(prot == i)
#    track2 = c(track2, track + 1)
#    prot = c(prot, i)
#  } else {
#    track2 = c(track2, 1)
#    prot = c(prot, i)
#  }
#}
#Data2$count = track2


new$pfam = as.character(new$pfam)


p = plot_ly(x = new$pfam, y = new$tpm, type = "scatter", mode = "markers", color = new$tpm, size = new$tpm, colors = "Blues")
p

p1 = plot_ly(x = new$pfam, y = new$tpm, type = "scatter", mode = "markers", color = new$newcolumn, size = new$tpm, colors = "Spectral")
p1


#save(new, p, file = "all3.RData")

load("all2.RData")

new2$pfam <- NULL
colSums(!is.na(new2))

new2[is.na(new2)] <- 0
colSums(new2)

new2 = new2[complete.cases(new2), ]



a = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
b = new2[is.na(new2$lepidophylla) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
c = new2[is.na(new2$wallacei) & is.na(new2$lepidophylla) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
d = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$lepidophylla) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
e = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$lepidophylla) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
f = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$lepidophylla) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
g = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$lepidophylla) & is.na(new2$stauntoniana),]
h = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$lepidophylla),]


trinity = c()
annot = c()

a5 = dim(a)[1]

for (i in seq(1:a5)) {
  a3 = a[i,1]
  a1 = grepl(a3, annotation$V5)
  a2 = which(a1)
  a4 = as.character(annotation[a2,1])
  a5 = as.character(annotation[a2,5])
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

write.table(trinity, file = "a1.txt", col.names = FALSE)
