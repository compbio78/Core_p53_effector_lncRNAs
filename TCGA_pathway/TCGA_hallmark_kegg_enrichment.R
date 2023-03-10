library('tidyr')
library('WebGestaltR')
library('ggplot2')
library(plyr)

##script 1: 
  ##Identifying lncRNA and cancer hallmark geneset associations across ten cancer types
  ##Identifying lncRNA and KEGG pathway associations across ten cancer types

regression_type<-"negative" #Uncomment for genes associated negatively with p53-effector lncRNAs
#regression_type<-"positive" #Uncomment for genes associated positively with p53-effector lncRNAs

cantype<-c('BLCA','BRCA','HNSC','KIRC','LIHC','LUAD','SKCM','STAD','OV','UCEC')
result<-lapply(cantype,function(c){
  df<-read.csv(paste0(c,'/lnc_mrna_',regression_type,'_regression.gmt'),check.names=FALSE,sep='\t',header=TRUE)
  res<-lapply(unique(df$lncrna),function(x){
    res<-NULL
    temp<-df[df$lncrna==x,]
    if(dim(temp)[1]>10){
      res<-WebGestaltR(enrichMethod="ORA",organism="hsapiens",
                       enrichDatabase="community-contributed_Hallmark50",
                       #enrichDatabase="pathway_KEGG", #uncomment for KEGG pathway analysis
                       enrichDatabaseType="genesymbol",interestGene=as.vector(temp$genesymbol),interestGeneType="genesymbol",
                       referenceGene=NULL,referenceSet="genome_protein-coding",referenceGeneType="genesymbol",
                       fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,reportNum=0,isoutput=FALSE,projectName="temp_delit/")
      res<-res[,c(1,2,8,9)]
      if(!is.null(res)){res$lncRNA_id<-as.character(x)}
      res
    } #end of if
  }) #end of inner lapply
  res1<-as.data.frame(do.call(rbind,res))
  res1$cantype<-c
  res1
}) #end of outer lapply
res1<-as.data.frame(do.call(rbind,result)) 
head(res1)
write.table(res1,file=paste0("hallmark_",regression_type,"_association_with_p53_effector_lncRNAs.csv"),sep="\t",quote=FALSE,row.names=FALSE) #replace hallmark with kegg for saving kegg pathway analysis results

###End of script1

##script 2: 
##Following script extracts source data for Fig. 4A
cantype<-c('BLCA','BRCA','HNSC','KIRC','LIHC','LUAD','SKCM','STAD','OV','UCEC')
df_neg<-read.csv("hallmark_negativeregress_res.csv",check.names=FALSE,sep='\t',header=TRUE) ## To have this file, first run script1 with hallmark data 
df_pos<-read.csv("hallmark_positiveregress_res.csv",check.names=FALSE,sep='\t',header=TRUE) ## To have this file, first run script1 with hallmark data

h_category<-read.csv("hallmark_categories.txt",sep="\t",header=TRUE,check.names=FALSE)
head(h_category)
h_category$Hallmark_Name<-paste('HALLMARK',h_category$Hallmark_Name,sep="_")
names(h_category)<-c('geneSet','category')

proc_dat<-function(df,category,corrtype){
  df$'comb'<-paste(df$lncRNA_id,df$geneSet,sep='.')
  df$corr<-corrtype
  df<-merge(df,category,by='geneSet',all=FALSE)
  df<-df[df$category=='proliferation',]
  result<-lapply(cantype,function(c){
    result1<-lapply(unique(df$geneSet),function(x){
      temp<-subset(df,df$geneSet==x & df$cantype ==c)
      c(as.character(x),as.character(c),length(unique(temp$lncRNA_id)))
    })
    result1<-as.data.frame(do.call(rbind,result1))
    result1$cantype<-c
    result1
  })
  df<-as.data.frame(do.call(rbind,result))
  names(df)<-c('geneSet','cantype','lnc_num','cantype1')
  df$comb<-paste(df$geneSet,df$cantype,sep='.')
  df
}

df1_neg<-proc_dat(df_neg,h_category,"neg")
df1_pos<-proc_dat(df_pos,h_category,"pos")
names(df1_neg)<-c('geneset_neg','cantype_neg','lnc_num_neg','cantype1_neg','comb')
names(df1_pos)<-c('geneset_pos','cantype_pos','lnc_num_pos','cantype1_pos','comb')
df1<-merge(df1_neg,df1_pos,by='comb',all=TRUE)
df1<-separate(df1,col=comb,into=c('geneset','cantype'),sep='\\.')
df1<-df1[,c('geneset','cantype','lnc_num_pos','lnc_num_neg')]
df1[is.na(df1)]<-0
write.table(df1,file='hallmark_lncnum_proliferation_suppress_induce_fig4A.csv',row.names=FALSE,quote=FALSE,sep=',')

##Following script extracts source data for Fig. 4B and Supp. Fig. S5
proc_dat<-function(df,category,corrtype){
  df$'comb'<-paste(df$lncRNA_id,df$geneSet,sep='.')
  df$corr<-corrtype
  df<-merge(df,category,by='geneSet',all=FALSE)
  df<-df[df$category=='proliferation',]
  df
}

df1_neg<-proc_dat(df_neg,h_category,"neg")
df1_pos<-proc_dat(df_pos,h_category,"pos")



proc_dat1<-function(df){
  res<-lapply(unique(df$lncRNA_id),function(x){
    temp<-subset(df,df$lncRNA_id==x)
    num_can<-length(unique(as.character(temp$cantype)))
    num_geneset<-length(unique(as.character(temp$geneSet)))
    c(as.character(x),num_can,num_geneset)
  })
  res1<-as.data.frame(do.call(rbind,res))
  res1
}


##proliferation suppressive lncRNAs
df1_neg_suppressive<-subset(df1_neg,df1_neg$geneSet !='HALLMARK_P53_PATHWAY')
df1_pos_suppressive<-subset(df1_pos,df1_pos$geneSet =='HALLMARK_P53_PATHWAY')
df1_suppressive<-rbind(df1_neg_suppressive,df1_pos_suppressive)
df1_suppressive<-proc_dat1(df1_suppressive)
names(df1_suppressive)<-c('lncrna_id','num_cantype_suppress_proliferation','num_suppressed_geneset')

##proliferation inducing lncRNAs
df1_neg_inducing<-subset(df1_neg,df1_neg$geneSet =='HALLMARK_P53_PATHWAY')
df1_pos_inducing<-subset(df1_pos,df1_pos$geneSet !='HALLMARK_P53_PATHWAY')
df1_inducing<-rbind(df1_neg_inducing,df1_pos_inducing)
df1_inducing<-proc_dat1(df1_inducing)
names(df1_inducing)<-c('lncrna_id','num_cantype_induce_proliferation','num_induced_geneset')

res<-merge(df1_suppressive,df1_inducing,by='lncrna_id',all=TRUE)
head(res)
res[is.na(res)]<-0
res$num_cantype_suppress_proliferation<-as.numeric(as.character(res$num_cantype_suppress_proliferation))
res$num_cantype_induce_proliferation<-as.numeric(as.character(res$num_cantype_induce_proliferation))
res$num_suppressed_geneset<-as.numeric(as.character(res$num_suppressed_geneset))
res$num_induced_geneset<-as.numeric(as.character(res$num_induced_geneset))

head(res)
res<-res[order(-res$num_cantype_suppress_proliferation,-res$num_suppressed_geneset),]
plot(seq(1,dim(res)[1]),res$num_cantype_suppress_proliferation,cex=res$num_suppressed_geneset/2,
     col=rgb(1,0,0,alpha=0.6),pch=16,ylim=c(1,10))
plot(seq(1,6),rep(4,6),cex=(c(0.5,1,1.5,2,2.5,3)),pch=16,col=rgb(1,0,0,alpha=0.6),xlim=c(0,7))
write.table(res,file='lncRNA_alter_proliferation_cancer_number_fig4b_supp_fig_S5.csv',row.names=FALSE,quote=FALSE,sep=',')

##End of script 2