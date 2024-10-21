library(ggplot2)
library(plotly)
library(plyr)

annotation = read.delim("ABIJ.tsv", header=F)
annotation2 = read.delim("JKAA.tsv", header=F)
annotation3 = read.delim("KJYC.tsv", header=F)
annotation4 = read.delim("KUXM.tsv", header=F)
annotation5 = read.delim("LGDQ.tsv", header=F)
annotation6 = read.delim("ZFGK.tsv", header=F)
annotation7 = read.delim("ZYCD.tsv", header=F)
annotation8 = read.delim("ZZOL.tsv", header=F)

#part A
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

#clean gene label
a<-gsub("(.*)_.*","\\1",Data$gene_id)
Data$gene_id <- a


justGeneP<-annotation[,c(1,5)]
justGeneTPM<-ABIJ.abundance[,c(1,5)]


colnames(Data)[colnames(Data)=="gene_id"] <- "target_id"
Data2 = merge(Data, justGeneTPM)
Data2 <- Data2[c(2,3)]
#just GO and TPM
Data2 = ddply(Data2, "go_id", numcolwise(sum))
#sum tpm values

Data3 = Data2

colnames(Data3)[colnames(Data3)=="tpm"] <- "lepidophylla"
Data2$newcolumn <- "lepidophylla"
#label, repeat for each sample - annotation2, annotation3


#Selaginella lepidophylla (ABIJ)
#Selaginella wallacei (JKAA)
#Selaginella willdenowii (KJYC)
#Selaginella selaginoides (KUXM)
#Selaginella apoda (LGDQ)
#Selaginella kraussiana (ZFGK)
#Selaginella acanthonota (ZYCD)
#Selaginella stauntoniana (ZZOL)


new2 = Data3
#use above on first run
#use below on subsequent runs
#new2 = merge.data.frame(new2, Data3, all = TRUE)


new = Data2
#use above on first run
#use below on subsequent runs
#new = rbind(new, Data2)


p0 = plot_ly(x = new$newcolumn, y = new$go_id, type = "scatter", mode = "markers", color = new$go_id, size = new$tpm, colors = "Spectral")

p0
#visualization


p1 = plot_ly(x = new$go_id, y = new$tpm, type = "scatter", mode = "markers", color = new$newcolumn, size = new$tpm, colors = "Spectral")

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
#unique GO terms

trinity = c()
annot = c()

a5 = dim(h)[1]
#change letter according to sample

for (i in seq(1:a5)) {
  a3 = h[i,1]
  a1 = grepl(a3, annotation8$V14)
  #change annotation variable according to sample
  a2 = which(a1)
  a4 = as.character(annotation8[a2,1])
  #change annotation variable according to sample
  a5 = as.character(annotation8[a2,14])
  #change annotation variable according to sample
  trinity = c(trinity, a4)
  annot = c(annot, a5)
}

annot = as.data.frame(annot)

trinity = as.data.frame(trinity)
clean <- gsub("(.*)_.*","\\1",trinity$trinity)
trinity$trinity <- clean
trinity2 = unique(trinity)
#fetch species-specific transcripts



#Load TSV file created by InterProScan as the "annotation" variable. Extract GO terms. Load TSV file "abundance" created by kallisto. Merge quantitative information with annotations. Create a bubble plot visualization. Fetch transcripts linked to species-specific GO terms.