##Load required R packages

library('tidyr')
library('WebGestaltR')
library('ggplot2')
library(plyr)
library(scatterplot3d)

##Script 1: Identifying lncRNA and cancer cell survival/growth associations across the cell lines of 10 cancer types

#Download 'Achilles_gene_effect.csv' file from depmap version 19Q3
#weblink: https://depmap.org/portal/download/all/

#Transpose the 'Achilles_gene_effect.csv' file and save it to 'transposed_Achilles_gene_effect.csv'
df<-read.csv('./crispr_19Q3/Achilles_gene_effect.csv',check.names=FALSE,sep=',',header=TRUE) 
df<-as.data.frame(t(df))
write.table(df,file=paste('./crispr_19Q3/transposed_Achilles_gene_effect.csv',sep=""),sep=',',quote=FALSE,row.names=TRUE,col.names=FALSE)

can_list<-c('BLCA','BRCA','HNSC','KIRC','LIHC','LUAD','OV','SKCM','STAD','UCEC')

df<-read.csv('./crispr_19Q3/transposed_Achilles_gene_effect.csv',check.names=FALSE,sep=',',header=TRUE)

regression_type<-"negative" #Uncomment for GSEA based on genesets associated negatively with p53-effector lncRNAs
#regression_type<-"positive" #Uncomment for GSEA based on genesets associated positively with p53-effector lncRNAs

lapply(can_list,function(cantype){
  dat<-df
  row.names(dat)<-dat[,1]
  names(dat)[1]<-'symbol'
  
  c<-read.csv(paste("./depmap/",cantype,"_mapped_cell_lines.txt",sep=""),check.names=FALSE,header=TRUE,sep="\t") #read list of cell lines of a  specific cancer type
  input<-dat[,names(dat) %in% c('symbol',as.character(rownames(c)))]
  input<-separate(input,col='symbol',into=c('symbol','id'),sep=" ")
  input<-input[,-2]
  #write.table(input,file=paste(cantype,"/",cantype,"_depmap_data.csv",sep=""),sep=',',quote=FALSE,row.names=FALSE)

  for (i in 2:length(input[1,])){
    x=input[,c(1,i)]
    x=subset(x,x[,2]!='NA')
    #GSEA using WebGestaltR
    res<-WebGestaltR(enrichMethod="GSEA",organism="others",
                     enrichDatabaseFile=paste0(cantype,"/lnc_mrna_",regression_type,"_regression.gmt"),
                     interestGene=x,fdrMethod="BH",sigMethod="fdr",fdrThr=1,isoutput=FALSE,
                     minNum=5, maxNum=10000,reportNum=0,perNum=1000,projectName=paste(cantype,regression_type,names(input[i]),sep="_"))
    
    res<-res[,-c(2,7,8,9,10)]
    res$cell_line<-names(input[i])
    if(i==2){
      write.table(res,file=paste(cantype,"/",regression_type,"_regression/",regression_type,"_regression_res_gsea_",cantype,".csv",sep=""),sep=',',quote=FALSE,row.names=FALSE,append=FALSE)
    }
    else{
      write.table(res,file=paste(cantype,"/",regression_type,"_regression/",regression_type,"_regression_res_gsea_",cantype,".csv",sep=""),sep=',',quote=FALSE,row.names=FALSE,append=TRUE)
    }
  }

}) #end of lapply

##script 1 end



##script 2: summary of all cancer types
fdr_cutoff<-0.001
meta_up_lnc_num<-49 #Number of predicted core p53-effector lncRNAs
cantype<-c('BLCA','BRCA','HNSC','KIRC','LIHC','LUAD','OV','SKCM','STAD','UCEC')
res_neg<-lapply(cantype,function(i){
  df<-read.csv(paste(i,"/negative_regression/negative_regression_res_gsea_",i,".csv",sep=''),
           check.names=FALSE,sep=',',header=TRUE,stringsAsFactors = FALSE)
  df$'cantype'<-i
  df<-subset(df,!df$geneSet=='geneSet')
  df
})

res_pos<-lapply(cantype,function(i){
  df<-read.csv(paste(i,"/positive_regression/positive_regression_res_gsea_",i,".csv",sep=''),
           check.names=FALSE,sep=',',header=TRUE,stringsAsFactors = FALSE)
  df$'cantype'<-i
  df<-subset(df,!df$geneSet=='geneSet')
  df
})

df_neg_reg<-as.data.frame(do.call("rbind",res_neg))
df_pos_reg<-as.data.frame(do.call("rbind",res_pos))

head(df_neg_reg)
tail(df_neg_reg)
dim(df_pos_reg)
dim(df_neg_reg)

df_pos_reg[,c(3,5)]<-lapply(df_pos_reg[,c(3,5)],as.character)
df_pos_reg[,c(3,5)]<-lapply(df_pos_reg[,c(3,5)],as.numeric)

df_neg_reg[,c(3,5)]<-lapply(df_neg_reg[,c(3,5)],as.character)
df_neg_reg[,c(3,5)]<-lapply(df_neg_reg[,c(3,5)],as.numeric)

df_pos_reg<-subset(df_pos_reg,df_pos_reg$normalizedEnrichmentScore < 0 & df_pos_reg$FDR < fdr_cutoff)
df_neg_reg<-subset(df_neg_reg,df_neg_reg$normalizedEnrichmentScore < 0 & df_neg_reg$FDR < fdr_cutoff)


dim(df_pos_reg)
dim(df_neg_reg)

df_pos_reg$'comb'<-paste(df_pos_reg$geneSet,df_pos_reg$cell_line,sep='.')
df_neg_reg$'comb'<-paste(df_neg_reg$geneSet,df_neg_reg$cell_line,sep='.')

df_pos_reg1<-subset(df_pos_reg,!(df_pos_reg$comb %in% df_neg_reg$comb))
df_neg_reg1<-subset(df_neg_reg,!(df_neg_reg$comb %in% df_pos_reg$comb))

id_file<-read.csv('cellline_id_to_name.txt',check.names=FALSE,sep='\t',header=TRUE)
df_pos_reg1<-merge(df_pos_reg1,id_file,by='cell_line',all = FALSE)
df_neg_reg1<-merge(df_neg_reg1,id_file,by='cell_line',all = FALSE)

posreg_dat<-df_pos_reg1
negreg_dat<-df_neg_reg1

write.table(posreg_dat,file=paste0("res_p53_effector_lnc_induce_growth_fdr_",fdr_cutoff,".csv"),sep='\t',quote=FALSE,row.names=FALSE,append=FALSE)
write.table(negreg_dat,file=paste0("res_p53_effector_lnc_suppress_growth_fdr_",fdr_cutoff,".csv"),sep='\t',quote=FALSE,row.names=FALSE,append=FALSE)

celllines<-unique(posreg_dat$cell_line)
res_pos<-lapply(celllines,function(x){
  temp<-posreg_dat[posreg_dat$cell_line==x,]
  c(as.character(x),dim(temp)[1])
})

res_pos<-as.data.frame(do.call(rbind,res_pos))
colnames(res_pos)<-c('celllines','tot_num_lnc_induce_growth')
#head(res_pos)
res_pos$tot_num_lnc_induce_growth<-as.numeric(res_pos$tot_num_lnc_induce_growth)

celllines<-unique(negreg_dat$cell_line)
res_neg<-lapply(celllines,function(x){
  temp<-negreg_dat[negreg_dat$cell_line==x,]
  c(as.character(x),dim(temp)[1])
})

res_neg<-as.data.frame(do.call(rbind,res_neg))
colnames(res_neg)<-c('celllines','tot_num_lnc_suppress_growth')
#head(res_neg)
res_neg$tot_num_lnc_suppress_growth<-as.numeric(res_neg$tot_num_lnc_suppress_growth)

merged_res<-merge(res_pos,res_neg,by='celllines',all=TRUE)
#dim(merged_res)
#head(merged_res)
merged_res$percent_lnc_induce_growth<-(merged_res$tot_num_lnc_induce_growth/meta_up_lnc_num)*100
merged_res$percent_lnc_suppress_growth <-(merged_res$tot_num_lnc_suppress_growth/meta_up_lnc_num)*100

head(merged_res)

temp_merged_res<-merged_res
temp_merged_res[is.na(temp_merged_res)]=0
wilcox_res<-wilcox.test(temp_merged_res$tot_num_lnc_induce_growth,temp_merged_res$tot_num_lnc_suppress_growth,alternative='two.sided')
wilcox_res$p.value
write.table(temp_merged_res,file=paste0("res_ts_onco_lnc_distribution_FDR_cutoff_",fdr_cutoff,".csv"),sep=',',row.names = TRUE,quote=FALSE)

#ECDF plots
ts<-temp_merged_res[,c('celllines','tot_num_lnc_suppress_growth')]
ts$category<-'ts'
names(ts)[2]<-'lnc_num'

onco<-temp_merged_res[,c('celllines',"tot_num_lnc_induce_growth")]
onco$category<-'onco'
names(onco)[2]<-'lnc_num'

ts_onco_lnc<-rbind(ts,onco)
p<- ggplot(ts_onco_lnc, aes(x=lnc_num, y=1-..y.., fill=category, colour=category)) +
  stat_ecdf(alpha=0.8,pad=TRUE,size=1)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank())

p + scale_colour_manual(values = c( 'grey2','red'))  


#3D scatter plots
plot_3d <- ddply(merged_res, .(merged_res$percent_lnc_suppress_growth, merged_res$percent_lnc_induce_growth), nrow)
plot_3d$num_celllines<-plot_3d$V1
head(plot_3d)
plot_3d[is.na(plot_3d)]=0
plot_3d=plot_3d[,-3]

scatterplot3d(plot_3d$num_celllines,
              plot_3d$`merged_res$percent_lnc_induce_growth`,
              plot_3d$`merged_res$percent_lnc_suppress_growth`,
              pch=16,color='palevioletred2',cex.axis=1,cex.symbol=1.2,angle=45,
              xlim=c(0,40),ylim=c(0,8),zlim=c(0,40),box=TRUE,grid=TRUE)

##script2 end

