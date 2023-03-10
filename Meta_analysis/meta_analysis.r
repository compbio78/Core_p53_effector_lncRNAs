##Following script identifies core p53-effector lncRNAs

library('RobustRankAggreg')
library('scatterplot3d')
gencode<-read.csv('gencode_v27.gtf',sep='\t',header=TRUE,check.names=FALSE)
gencode<-subset(gencode,gencode$lnc_ensid!='ENSG00000198788') #Removed MUC2 protein-coding gene 

chipseq_p53_tar<-read.csv('p53-target_lnc_chipbase_atleast_5_data.txt',header=TRUE,sep=',',check.names=FALSE)
rownames(chipseq_p53_tar)<-chipseq_p53_tar$lncRNA

path="input_RNAseq/"
filenames <- list.files(path)  
filenames

fc <- lapply(filenames,function(x){
  data_name<-stringr::str_sub(x, end=-5)
  data<-read.csv(paste0(path,x), header=TRUE, check.names=FALSE,sep=',')
  names(data)[1]=paste(data_name,"logFC",sep='_')
  data
})

merged_fc<-Reduce(function(x, y) merge(x, y, by = "ensid", all = TRUE), fc)
row.names(merged_fc)<-merged_fc$ensid
merged_fc<-merged_fc[,-1]
head(merged_fc)

na_rows <- apply(merged_fc, 1, function(x) all(is.na(x)))
merged_fc<-merged_fc[ !na_rows, ]
merged_fc$upregulation_count<-rowSums(merged_fc>0,na.rm=TRUE)
merged_fc$downregulation_count<-rowSums(merged_fc<0,na.rm=TRUE)
merged_fc$upregulation_proportion<-merged_fc$upregulation_count/(merged_fc$upregulation_count+merged_fc$downregulation_count)
merged_fc$downregulation_proportion<-merged_fc$downregulation_count/(merged_fc$upregulation_count+merged_fc$downregulation_count)
consistent_up<-merged_fc[merged_fc$upregulation_proportion > 0.667,]
consistent_down<-merged_fc[merged_fc$downregulation_proportion > 0.667,]


#RobustRankAggreg
rra<-function(filenames,type){
  mylist <- lapply(filenames,function(x){
    data_name<-stringr::str_sub(x, end=-5)
    data<-read.csv(paste0(path,x), header=TRUE, check.names=FALSE,sep=',')
    if(type=='up'){
      data<-data[order(data$logFC,decreasing=TRUE),]
      data<-rownames(data[data$logFC > 0,])
    }
    
    if(type=='down'){
      data<-data[order(data$logFC),]
      data<-rownames(data[data$logFC < 0,])
    }
    data
  })
  
  tot_genes<-length(unique(unlist(mylist, recursive = FALSE)))
  
  r = rankMatrix(mylist, N = tot_genes, full=TRUE)
  ar = aggregateRanks(rmat = r)
  
  ar_cv<-lapply(c(1:ncol(r)),function(x){
    aggregateRanks(rmat = r[,-x])
  })
  
  merged_ar_cv<-Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), ar_cv)
  merged_ar<-merge(ar,merged_ar_cv,by="Name",all=TRUE)
  row.names(merged_ar)<-merged_ar[,1]
  merged_ar<-merged_ar[,-1]
  merged_ar$p_val<-apply(merged_ar,1,mean,na.rm=TRUE)
  merged_ar$fdr<-p.adjust(as.numeric(as.character(merged_ar$p_val)),method="BH",n=length(merged_ar$p_val))
  
  if(type=='up')
    {
    merged_ar<-subset(merged_ar,round(merged_ar$fdr,3)<=0.050 & rownames(merged_ar) %in% rownames(consistent_up) & 
                        rownames(merged_ar) %in% gencode$lnc_ensid & rownames(merged_ar) %in% rownames(chipseq_p53_tar))
    }
  if(type=='down')
    {
    merged_ar<-subset(merged_ar,round(merged_ar$fdr,3)<=0.050 & rownames(merged_ar) %in% rownames(consistent_down) & 
                        rownames(merged_ar) %in% gencode$lnc_ensid & rownames(merged_ar) %in% rownames(chipseq_p53_tar))
    }
   
  merged_ar<-merge(merged_ar,chipseq_p53_tar,by='row.names',all=FALSE)
  rownames(merged_ar)<-merged_ar$Row.names
  merged_ar<-merged_ar[,c('fdr','num_chipseq')]  
  merged_ar
}

p53_lnc_up<-rra(filenames,'up')
p53_lnc_down<-rra(filenames,'down')

##Merge fold-change data
p53_lnc_up_fc<-merge(p53_lnc_up,merged_fc[,seq(1,14)],by='row.names')
rownames(p53_lnc_up_fc)<-p53_lnc_up_fc$Row.names
p53_lnc_up_fc<-p53_lnc_up_fc[,-1]
p53_lnc_up_fc<-p53_lnc_up_fc[order(p53_lnc_up_fc$fdr),]
write.table(p53_lnc_up_fc,file="res_rra.csv",row.names = TRUE, quote=FALSE, sep=',')

##3d scatter plot - Fig. 2E
dat<-merge(p53_lnc_up,merged_fc[,seq(1,15)],by='row.names')
rownames(dat)<-dat$Row.names
dat<-dat[,-1]
head(dat)
dat$median_fc<-apply(dat,1,function(x){
  median(x[3:16],na.rm=TRUE)
})
dat$fdr_log10<-(-1)*log(dat$fdr,10)
dat1<-dat[rownames(dat) %in% c("ENSG00000224294","ENSG00000225511","ENSG00000250337","ENSG00000234546","ENSG00000248429","ENSG00000253878","ENSG00000246640","ENSG00000182165","ENSG00000263718"),] #known p53-effector lncRNAs
dat2<-dat[!(rownames(dat) %in% c("ENSG00000224294","ENSG00000225511","ENSG00000250337","ENSG00000234546","ENSG00000248429","ENSG00000253878","ENSG00000246640","ENSG00000182165","ENSG00000263718")),]  #Remaining lncRNAs


dat1$col<-rgb(.01,.9,.1)
dat2$col<-rgb(0.9,.01,.3,alpha=0.4)
dat<-rbind(dat1,dat2)
   
scatterplot3d(dat$median_fc,dat$upregulation_count,dat$fdr_log10,pch=16,
              color=dat$col,cex.axis=1.1,cex.symbol=1.5,angle=30,label.tick.marks = T)


##End of script
