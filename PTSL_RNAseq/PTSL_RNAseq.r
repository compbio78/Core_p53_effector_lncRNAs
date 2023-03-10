library(edgeR)
library(limma)
library(grDevices)
library("org.Hs.eg.db")
library("ggplot2")
library('tidyr')
library(data.table)

##Script 1
## Differentially expressed gene identification using edgeR

#RNA-sequencing was performed to evaluate gene expression changes following ectopic expression of PTSL or vector control in Lung Adenocarcinoma cells.  
#The RNA-seq data can be downloaed from Gene Expression Omnibus (Accession ID is provided in the article).

#prepare following raw read counts table and save to the file "PTSL_rnaseq_read_count.txt"
  #Ensembl_Gene_ID control_rep1 control_rep2 control_rep3 control_rep4 ptsl_rep1 ptsl_rep2 ptsl_rep3 ptsl_rep4
#1 ENSG00000282222            0            0            0            0         0         0         0         0
#2 ENSG00000282221           39           26           32           22        16        13        12         0
#3 ENSG00000212040            0            0            0            0         0         0         0         0
#4 ENSG00000110514         4019         3168         3004         3668      3786      3102      5948      1298
#5 ENSG00000086015         5344         5119         4545         4803      4265      3193      6688      1426
#6 ENSG00000207355            0            0            0            0         0         0         0         0

d <-read.csv(file="PTSL_rnaseq_read_count.txt", header=TRUE,sep="\t",check.names=FALSE)
head(d)
rownames(d)<-d[,1]
d<-d[,-1]

colnames( d ) <- paste(c(rep("ctrl",4),rep("expt",4)),c(1:4,1:4),sep="")
colnames(d)

group <- c(rep("ctrl",4) , rep("expt",4))
print(group)

y <- DGEList(counts=d, group = group)
dim(y)

keep1<-rowSums(cpm(y[,1:8]) >= 1)>=4
y<-y[keep1,]
dim(y)
y <- calcNormFactors(y)

y <- estimateCommonDisp(y, verbose=TRUE)
y<-estimateTagwiseDisp(y)
et <- exactTest(y,pair=c("ctrl","expt"))
topTags(et,n=5)
summary(decideTestsDGE(et, adjust.method="BH"))

resultsTbl.tgw <- topTags(et,n=nrow(et$table))$table
head(resultsTbl.tgw)
write.table(resultsTbl.tgw, file="1_res_de.csv", sep = "," , row.names = TRUE)

resultsTbl.tgw$symbol<-mapIds(org.Hs.eg.db,keys=row.names(resultsTbl.tgw),column="SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(resultsTbl.tgw, file="2_symbol_res.csv", sep = "," , row.names = TRUE,quote=F)

no_na_resultsTbl.tgw<-resultsTbl.tgw[!is.na(resultsTbl.tgw$symbol),]
write.table(no_na_resultsTbl.tgw, file="3_symbol_NA_removed_res.csv", sep = "," , row.names = TRUE,quote=F)

head(no_na_resultsTbl.tgw)
gsea_input<-no_na_resultsTbl.tgw[,c('symbol','logFC')]
write.table(gsea_input, file="4_gsea_input.rnk", sep = "\t" , row.names = FALSE,quote=F)


##Script 2
##code for Figure 5c
data <-read.csv(file="3_symbol_NA_removed_res.csv", header=TRUE,sep=",",check.names=FALSE)  ##Get this file from Script 1
head(data)

#get TCGA Lung Adenocarcinoma lncRNA-mRNA expression association data
negative<-read.csv('regression_tcga_luad/lnc_mrna_negative_regression.gmt',header=TRUE,sep="\t",check.names=FALSE)
head(negative)
negative<-negative[negative$geneset=='ENSG00000253878',] #geneset variable holds lncRNA ensembl id; #ENSG00000253878=PTSL
negative_fc<-subset(data,data$symbol %in% negative$genesymbol)
head(negative_fc)
dim(negative_fc)

positive<-read.csv('regression_tcga_luad/lnc_mrna_positive_regression.gmt',header=TRUE,sep="\t",check.names=FALSE)
head(positive)
positive<-positive[positive$geneset=='ENSG00000253878',] #ENSG00000253878=PTSL
dim(positive)
positive_fc<-subset(data,data$symbol %in% positive$genesymbol)
head(positive_fc)
dim(positive_fc)

plot(ecdf(as.numeric(as.character(negative_fc$logFC))),verticals= TRUE,xlim=c(-1,1),ylim=c(0,1),cex=0.5,axes=T,
     lty=1, pch=19, col="blue", lwd= 2, main=" ", xlab="",ylab="")

par(new=T)

plot(ecdf(as.numeric(as.character(positive_fc$logFC))),verticals= TRUE,xlim=c(-1,1),ylim=c(0,1),cex=0.5,axes=T,
     lty=1, pch=19, col="red", lwd= 2, main=" ", xlab="",ylab="")

p<-wilcox.test(as.numeric(as.character(negative_fc$logFC)),as.numeric(as.character(positive_fc$logFC)), alternative='two.sided')
p$p.value

##Script 3
##code for Figure 5D
data <-read.csv(file="3_symbol_NA_removed_res.csv", header=TRUE,sep=",",check.names=FALSE)  ##Get this file from Script 1
head(data)

tsgene<-read.csv('ts_genes.txt')
oncogene<-read.csv('onco_genes.txt')

tsgene_data<-subset(data,data$symbol %in% tsgene$gene_symbol)
oncogene_data<-subset(data,data$symbol %in% oncogene$gene_symbol)
oncogene_data<-oncogene_data[!(rownames(oncogene_data)=='ENSG00000274512'),] #ENSG00000274512=TBC1D3L, which was incorrectly annotated as TBC1D3 in org.Hs.eg.db; hence excluded


plot(ecdf(as.numeric(as.character(oncogene_data$logFC))),verticals= TRUE,ylim=c(0,1),xlim=c(-1,1),cex=0.4,axes=T,
     lty=1, pch=19, col="purple", lwd= 3, main=" ", xlab="",ylab="")

par(new=T)

plot(ecdf(as.numeric(as.character(tsgene_data$logFC))),verticals= TRUE,ylim=c(0,1),xlim=c(-1,1),cex=0.4,axes=T,
     lty=1, pch=19, col="green", lwd= 3, main=" ", xlab="",ylab="")


p<-wilcox.test(as.numeric(as.character(oncogene_data$logFC)),as.numeric(as.character(tsgene_data$logFC)), alternative='two.sided')
p$p.value

##script 4
#code for Figure 5E
data <-read.csv(file="3_symbol_NA_removed_res.csv", header=TRUE,sep=",",check.names=FALSE)  ##Get this file from Script 1
head(data)

favorable<-read.csv("tcga_luad_favorable_genes.txt")
unfavorable<-read.csv("tcga_luad_unfavorable_genes.txt")

favorable_data<-subset(data,data$symbol %in% favorable$gene_symbol)
unfavorable_data<-subset(data,data$symbol %in% unfavorable$gene_symbol)

plot(ecdf(as.numeric(as.character(favorable_data$logFC))),verticals= TRUE,ylim=c(0,1),xlim=c(-1,1),cex=0.3,axes=T,
     lty=1, pch=19, col="orange", lwd= 2, main=" ", xlab="",ylab="")

par(new=T)

plot(ecdf(as.numeric(as.character(unfavorable_data$logFC))),verticals= TRUE,ylim=c(0,1),xlim=c(-1,1),cex=0.5,axes=T,
     lty=1, pch=19, col="black", lwd= 3, main=" ", xlab="",ylab="")


p<-wilcox.test(as.numeric(as.character(favorable_data$logFC)),as.numeric(as.character(unfavorable_data$logFC)), alternative='two.sided')
p$p.value

##End


