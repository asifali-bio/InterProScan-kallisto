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


#part A
justGeneP<-annotation8[,c(1,5)]
justGeneTPM<-ZZOL.abundance[,c(1,5)]

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

colnames(Data3)[colnames(Data3)=="tpm"] <- "stauntoniana"
Data2$newcolumn <- "stauntoniana"
#label


#Selaginella lepidophylla (ABIJ)
#Selaginella wallacei (JKAA)
#Selaginella willdenowii (KJYC)
#Selaginella selaginoides (KUXM)
#Selaginella apoda (LGDQ)
#Selaginella kraussiana (ZFGK)
#Selaginella acanthonota (ZYCD)
#Selaginella stauntoniana (ZZOL)
#new2 = Data3
new2 = merge.data.frame(new2, Data3, all = TRUE)


#new = Data2
new = rbind(new, Data2)
#combine all data


p = plot_ly(x = new$pfam, y = new$tpm, type = "scatter", mode = "markers", color = new$newcolumn, size = new$tpm, colors = "Spectral")

p


p1 = plot_ly(x = new$newcolumn, y = new$pfam, type = "scatter", mode = "markers", color = new$pfam, size = new$tpm, colors = "Spectral")

p1



#part B
a = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
b = new2[is.na(new2$lepidophylla) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
c = new2[is.na(new2$wallacei) & is.na(new2$lepidophylla) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
d = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$lepidophylla) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
e = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$lepidophylla) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
f = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$lepidophylla) & is.na(new2$acanthonota) & is.na(new2$stauntoniana),]
g = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$lepidophylla) & is.na(new2$stauntoniana),]
h = new2[is.na(new2$wallacei) & is.na(new2$willdenowii) & is.na(new2$selaginoides) & is.na(new2$apoda) & is.na(new2$kraussiana) & is.na(new2$acanthonota) & is.na(new2$lepidophylla),]
#unique protein domains


annotation <- annotation[which(annotation$V9 < 0.001),]
annotation2 <- annotation2[which(annotation2$V9 < 0.001),]
annotation3 <- annotation3[which(annotation3$V9 < 0.001),]
annotation4 <- annotation4[which(annotation4$V9 < 0.001),]
annotation5 <- annotation5[which(annotation5$V9 < 0.001),]
annotation6 <- annotation6[which(annotation6$V9 < 0.001),]
annotation7 <- annotation7[which(annotation7$V9 < 0.001),]
annotation8 <- annotation8[which(annotation8$V9 < 0.001),]
#filter by e-value


trinity = c()
annot = c()

a5 = dim(a)[1]
#change letter according to sample

for (i in seq(1:a5)) {
  a3 = a[i,1]
  a1 = grepl(a3, annotation$V5)
  #change annotation variable according to sample
  a2 = which(a1)
  a4 = as.character(annotation[a2,1])
  #change annotation variable according to sample
  a5 = as.character(annotation[a2,5])
  #change annotation variable according to sample
  trinity = c(trinity, a4)
  annot = c(annot, a5)
}

annot = as.data.frame(annot)

trinity = as.data.frame(trinity)
clean <- gsub("(.*)_.*","\\1",trinity$trinity)
trinity$trinity <- clean
trinity2 = unique(trinity)



#Similar to the previous step, we merge the protein domain annotations by InterProScan via Pfam with quantitative information at the gene-level by kallisto. We plot the data. We then fetch the transcripts associated with the species-specific Pfam protein domains filtered by E-value.