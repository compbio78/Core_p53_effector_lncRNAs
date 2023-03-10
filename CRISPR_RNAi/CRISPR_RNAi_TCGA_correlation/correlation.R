library('tidyr')
library('WebGestaltR')
library('ggplot2')
library(plyr)

#Script1 - Correlation between CRISPR and RNAi (Fig. 2G)

crispr<-read.csv("crispr_res_p53_effector_lnc_suppress_growth_fdr_0.001.csv",check.names=FALSE,sep='\t',header=TRUE) #Generated this file from CRISPR data analysis
rnai<-read.csv("rnai_res_p53_effector_lnc_suppress_growth_fdr_0.001.csv",check.names=FALSE,sep='\t',header=TRUE) ##Generated this file from RNAi data analysis

crispr_freq<-as.data.frame(table(crispr$geneSet))
names(crispr_freq)[1]<-'lncid'
crispr_freq$percent<-(crispr_freq$Freq/255)*100 #total 255 celllines were taken for the analysis

rnai_freq<-as.data.frame(table(rnai$geneSet))
names(rnai_freq)[1]<-'lncid'
rnai_freq$percent<-(rnai_freq$Freq/311)*100 #total 311 celllines were taken for the analysis

merged_dat<-merge(crispr_freq,rnai_freq,by='lncid',all=TRUE)
names(merged_dat)<-c('lncid','crispr_num_cell_lines','crispr_percent_cell_lines','rnai_num_cell_lines','rnai_percent_cell_lines')

x<-cor.test(merged_dat$crispr_percent_cell_lines,merged_dat$rnai_percent_cell_lines,method='spearman')
x$p.value
x$estimate


ggplot(merged_dat, aes(x=crispr_percent_cell_lines,y=rnai_percent_cell_lines)) +
  geom_point(size=2,fill=rgb(0.2,0.8,0.2),shape=16,colour=3)+ 
  geom_smooth(method=lm) +
  xlim(0,80) + ylim(0,80) +
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))

merged_dat<-merged_dat[order(merged_dat$crispr_percent_cell_lines,decreasing=TRUE),]

##Script2 - Correlation between crispr rnai and hallmark (Fig. 2H)

crispr_rnai<-merged_dat
rownames(crispr_rnai)<-crispr_rnai$lncid

hallmark<-read.csv("lncRNA_alter_proliferation_cancer_number.csv",sep=',',header=TRUE) 
head(hallmark)
hallmark<-hallmark[,c('lncrna_id','num_cantype_suppress_proliferation')]
hallmark$'num_cantype_suppress_proliferation_percent'<-(hallmark$num_cantype_suppress_proliferation/10)*100
rownames(hallmark)<-hallmark$lncrna_id


df<-merge(hallmark,crispr_rnai,by="row.names",all=TRUE)
head(df)
rownames(df)<-df$Row.names

##########################measure correlation between crispr and hallmark and between RNAi and hallmark
#RNAi-Hallmark
x<-cor.test(df$num_cantype_suppress_proliferation_percent,df$rnai_percent_cell_lines,method='spearman',na.rm=TRUE)
x$p.value
x$estimate

#crispr-hallmark
x<-cor.test(df$num_cantype_suppress_proliferation_percent,df$crispr_percent_cell_lines,method='spearman',na.rm=TRUE)
x$p.value
x$estimate

df$rnai_percent_cell_lines<- (-1)*df$rnai_percent_cell_lines

rnai<-df[,c("rnai_percent_cell_lines","num_cantype_suppress_proliferation_percent")]
rnai$type<-'rnai'
names(rnai)[1]<-'rnai_crispr'
rnai$color<-"rgb(0.1,0.1,0.9)"

crispr<-df[,c("crispr_percent_cell_lines","num_cantype_suppress_proliferation_percent")]
crispr$type<-'crispr'
names(crispr)[1]<-'rnai_crispr'
crispr$color<-"rgb(0.9,0.1,0.1)"


df1<-rbind(rnai,crispr)
head(df1)


ggplot(df1, aes(x=rnai_crispr,y=num_cantype_suppress_proliferation_percent)) +
  geom_point(size=2,fill=rgb(0.1,0.8,0.1),shape=16,colour=3,alpha=1)+ 
  geom_smooth(method=loess) +
  xlim(-80,80) + ylim(0,100) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))

##End of script 1