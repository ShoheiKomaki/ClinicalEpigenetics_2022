# characterize CpGs

# load packages
library(reshape2)
library(ggplot2)
library(data.table)
library(coin)
library(missMethyl)

# load datasets
df1=read.table("path/to/SDvsCoef.tsv",header=T,stringsAsFactors=F,sep="\t")
tmp1=as.data.frame(fread("path/to/HM450manifestfile"),stringsAsFactors=F)
df2=tmp1[,c("probeID","designType")]

# merge
tmp1=merge(df1,df2,by.x="cgID",by.y="probeID",all.x=T,all.y=F)
df4=tmp1

# calc contribution to the DNAm age fluctuation
df4$contribution=df4$coef * df4$meanSD
df4$absContribution=abs(df4$contribution)
g=ggplot(df4,aes(x=designType,y=log10(absContribution)))+geom_violin()+geom_boxplot(outlier.shape=NA,width=0.2)+facet_wrap(~clock,scales="free")+theme_classic()+theme(strip.background=element_blank(),strip.text.x=element_text(size=15))+xlab("Probe type")+ylab(expression(paste(log[10],"(Contribution to fluctuation)")))
ggsave("contribution2DNAmAgeFluctuation.pdf",g,width=8,height=3)

# compare contribution of probe types
tmp1=as.data.frame(ftable(df4[,c("clock","designType")]))
tmp2=tmp1[order(tmp1$clock),]
tmp2
tmp3=NULL
for(clock in c("Horvath 2013","Horvath 2018","Levine 2018")){
 tmpA=df4[df4$clock==clock,]
 tmpB=tmpA[tmpA$designType=="I","absContribution"]
 tmpC=tmpA[tmpA$designType=="II","absContribution"]
 tmpD=length(tmpB)
 tmpE=length(tmpC)
 tmpF=wilcox_test(absContribution~factor(designType),data=tmpA)
 tmpG=pvalue(tmpF)
 tmp3=rbind(tmp3,c(clock,tmpD,tmpE,mean(tmpB),mean(tmpC),tmpG))
}
tmp4=as.data.frame(tmp3,stringsAsFactors=F)
colnames(tmp4)=c("clock","nTypeI","nTypeII","meanTypeI","meanTypeII","wilcoxP")
write.table(tmp4,"probeTypeComparison.tsv",row.names=F,quote=F,sep="\t")
df5=tmp4 
}

# enrichment analyses
tmp1=unique(df3[,c("probe_ID","trait")])
for(clock in c("Horvath 2013","Horvath 2018","Levine 2018")){
 print(clock)
 tmpA=df4[df4$clock==clock,]
 tmpA$absContribution=abs(tmpA$contribution)
 tmpB=tmpA[tmpA$absContribution>quantile(tmpA$absContribution,prob=0.75),"cgID"]
 tmpC=tmpA[tmpA$absContribution<quantile(tmpA$absContribution,prob=0.25),"cgID"]
 tmpHighGO=gometh(sig.cpg=tmpB,all.cpg=tmpA$cgID,collection="GO")
 tmpHighGO=tmpHighGO[tmpHighGO$FDR<=0.05,]
 print(paste0("High contribution: GO: ",nrow(tmpHighGO)))
 tmpLowGO=gometh(sig.cpg=tmpC,all.cpg=tmpA$cgID,collection="GO")
 tmpLowGO=tmpLowGO[tmpLowGO$FDR<=0.05,]
 print(paste0("Low contribution: GO: ",nrow(tmpLowGO)))
 tmpHighKEGG=gometh(sig.cpg=tmpB,all.cpg=tmpA$cgID,collection="KEGG")
 tmpHighKEGG=tmpHighKEGG[tmpHighKEGG$FDR<=0.05,]
 print(paste0("High contribution: KEGG: ",nrow(tmpHighKEGG)))
 tmpLowKEGG=gometh(sig.cpg=tmpB,all.cpg=tmpA$cgID,collection="KEGG")
 tmpLowKEGG=tmpLowKEGG[tmpLowKEGG$FDR<=0.05,]
 print(paste0("Low contribution: KEGG: ",nrow(tmpLowKEGG)))
}
# nothing enriched

# regional annotations (group by SD)
tmp1=as.data.frame(fread("path/to/HM450/annotationFile"),stringsAsFactors=F)
df6=tmp1[,c("cgID","annot.type")]
tmp2=data.frame(annot=unique(df6$"annot.type"),stringsAsFactors=F)
for(clock in c("Horvath 2013","Horvath 2018","Levine 2018")){
 tmpA=df4[df4$clock==clock,]
 tmpB=tmpA[tmpA$meanSD>quantile(tmpA$meanSD,prob=0.75),"cgID"]
 tmpC=tmpA[tmpA$meanSD<quantile(tmpA$meanSD,prob=0.25),"cgID"]
 tmpD=as.data.frame(ftable(df6[df6$cgID%in%tmpB,"annot.type"]),stringsAsFactors=F)
 tmpE=as.data.frame(ftable(df6[df6$cgID%in%tmpC,"annot.type"]),stringsAsFactors=F)
 tmpF=merge(tmpD,tmpE,by="Var1",all=T)
 tmpG=as.data.frame(ftable(df6[df6$cgID%in%tmpA$cgID,"annot.type"]),stringsAsFactors=F)
 tmpH=merge(tmpG,tmpF,by="Var1",all=T)
 colnames(tmpH)=c("annot",paste0(clock,"_all"),paste0(clock,"_var"),paste0(clock,"_stbl"))
 tmp2=merge(tmp2,tmpH,by="annot",all=T)
 ## Pearson's Chi-squared test
 a=tmpH[grep("cpg",tmpH$annot),c(1,3,4)]
 rownames(a)=a$annot
 b=na.omit(a[,-1])
 c=as.numeric(apply(b,1,sum))
 d=b[c>10,]
 print(paste0(clock,": CpG annotations"))
 print(chisq.test(d)$p.value)
 a=tmpH[grep("gene",tmpH$annot),c(1,3,4)]
 rownames(a)=a$annot
 b=na.omit(a[,-1])
 c=as.numeric(apply(b,1,sum))
 d=b[c>10,]
 print(paste0(clock,": Genic annotations"))
 print(chisq.test(d)$p.value)
}

## cpg annotations
tmp3=tmp2[grep("cpg",tmp2$annot),]
for(c in 2:ncol(tmp3)){
 tmpA=as.numeric(tmp3[,c])
 tmpA[is.na(tmpA)]=0
 tmpB=100*tmpA/sum(tmpA)
 tmp3[,c]=tmpB
}
tmp4=melt(tmp3,id.vars=c("annot"))
tmp4$variable=as.character(tmp4$variable)
tmp4$annot=sub("inter","open-sea",sub("hg19_cpg_","",tmp4$annot))
tmp4$clock=sub("_.*$","",tmp4$variable)
tmp4$category=sub("^.*_","",tmp4$variable)
tmp4$annot=factor(tmp4$annot,levels=c("open-sea","shelves","shores","islands"))
g=ggplot(tmp4,aes(x=category,y=value,fill=annot))+geom_bar(position="stack",stat="identity",colour="white")+coord_flip()+theme_bw()+ylab("Frequency (%)")+xlab("")+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=10))+theme_classic()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),strip.background=element_blank(),strip.text.x=element_text(size=15))+facet_wrap(~clock)
ggsave("CpGAnnotationsBySD+legend_25percent.pdf",g,width=8,height=8)
g=g+theme(legend.position="NONE")
ggsave("CpGAnnotationsBySD_25percent.pdf",g,width=8,height=1.5)
## genic annotations
tmp3=tmp2[grep("cpg",tmp2$annot,invert=T),]
for(c in 2:ncol(tmp3)){
 tmpA=as.numeric(tmp3[,c])
 tmpA[is.na(tmpA)]=0
 tmpB=100*tmpA/sum(tmpA)
 tmp3[,c]=tmpB
}
tmp4=melt(tmp3,id.vars=c("annot"))
tmp4$variable=as.character(tmp4$variable)
tmp4$annot=sub("hg19_genes_","",tmp4$annot)
tmp4$clock=sub("_.*$","",tmp4$variable)
tmp4$category=sub("^.*_","",tmp4$variable)
tmp4$annot=factor(tmp4$annot,levels=rev(c("exons","introns","intronexonboundaries","5UTRs","3UTRs","promoters","1to5kb","intergenic")))
g=ggplot(tmp4,aes(x=category,y=value,fill=annot))+geom_bar(position="stack",stat="identity",colour="white")+coord_flip()+theme_bw()+ylab("Frequency (%)")+xlab("")+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=10))+theme_classic()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),strip.background=element_blank(),strip.text.x=element_text(size=15))+facet_wrap(~clock)
ggsave("GenicAnnotationsBySD+legend_25percent.pdf",g,width=8,height=8)
g=g+theme(legend.position="NONE")
ggsave("GenicAnnotationsBySD_25percent.pdf",g,width=8,height=1.5)


# regional annotations (group by contribution)
tmp1=as.data.frame(fread("path/to/HM450/annotationFile"),stringsAsFactors=F)
df6=tmp1[,c("cgID","annot.type")]
tmp2=data.frame(annot=unique(df6$"annot.type"),stringsAsFactors=F)
for(clock in c("Horvath 2013","Horvath 2018","Levine 2018")){
 tmpA=df4[df4$clock==clock,]
 tmpA$absContribution=abs(tmpA$contribution)
 tmpB=tmpA[tmpA$absContribution>quantile(tmpA$absContribution,prob=0.75),"cgID"]
 tmpC=tmpA[tmpA$absContribution<quantile(tmpA$absContribution,prob=0.25),"cgID"]
 tmpD=as.data.frame(ftable(df6[df6$cgID%in%tmpB,"annot.type"]),stringsAsFactors=F)
 tmpE=as.data.frame(ftable(df6[df6$cgID%in%tmpC,"annot.type"]),stringsAsFactors=F)
 tmpF=merge(tmpD,tmpE,by="Var1",all=T)
 tmpG=as.data.frame(ftable(df6[df6$cgID%in%tmpA$cgID,"annot.type"]),stringsAsFactors=F)
 tmpH=merge(tmpG,tmpF,by="Var1",all=T)
 colnames(tmpH)=c("annot",paste0(clock,"_all"),paste0(clock,"_high"),paste0(clock,"_low"))
 tmp2=merge(tmp2,tmpH,by="annot",all=T)
 ## Pearson's Chi-squared test
 a=tmpH[grep("cpg",tmpH$annot),c(1,3,4)]
 rownames(a)=a$annot
 b=na.omit(a[,-1])
 c=as.numeric(apply(b,1,sum))
 d=b[c>10,]
 print(paste0(clock,": CpG annotations"))
 print(chisq.test(d)$p.value)
 a=tmpH[grep("gene",tmpH$annot),c(1,3,4)]
 rownames(a)=a$annot
 b=na.omit(a[,-1])
 c=as.numeric(apply(b,1,sum))
 d=b[c>10,]
 print(paste0(clock,": Genic annotations"))
 print(chisq.test(d)$p.value)
}

## cpg annotations
tmp3=tmp2[grep("cpg",tmp2$annot),]
for(c in 2:ncol(tmp3)){
 tmpA=as.numeric(tmp3[,c])
 tmpA[is.na(tmpA)]=0
 tmpB=100*tmpA/sum(tmpA)
 tmp3[,c]=tmpB
}
tmp4=melt(tmp3,id.vars=c("annot"))
tmp4$variable=as.character(tmp4$variable)
tmp4$annot=sub("inter","open-sea",sub("hg19_cpg_","",tmp4$annot))
tmp4$clock=sub("_.*$","",tmp4$variable)
tmp4$category=sub("^.*_","",tmp4$variable)
tmp4$annot=factor(tmp4$annot,levels=c("open-sea","shelves","shores","islands"))
g=ggplot(tmp4,aes(x=category,y=value,fill=annot))+geom_bar(position="stack",stat="identity",colour="white")+coord_flip()+theme_bw()+ylab("Frequency (%)")+xlab("")+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=10))+theme_classic()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),strip.background=element_blank(),strip.text.x=element_text(size=15))+facet_wrap(~clock)
ggsave("CpGAnnotationsByContribution+legend_25percent.pdf",g,width=8,height=8)
g=g+theme(legend.position="NONE")
ggsave("CpGAnnotationsByContribution_25percent.pdf",g,width=8,height=1.5)
## genic annotations
tmp3=tmp2[grep("cpg",tmp2$annot,invert=T),]
for(c in 2:ncol(tmp3)){
 tmpA=as.numeric(tmp3[,c])
 tmpA[is.na(tmpA)]=0
 tmpB=100*tmpA/sum(tmpA)
 tmp3[,c]=tmpB
}
tmp4=melt(tmp3,id.vars=c("annot"))
tmp4$variable=as.character(tmp4$variable)
tmp4$annot=sub("hg19_genes_","",tmp4$annot)
tmp4$clock=sub("_.*$","",tmp4$variable)
tmp4$category=sub("^.*_","",tmp4$variable)
tmp4$annot=factor(tmp4$annot,levels=rev(c("exons","introns","intronexonboundaries","5UTRs","3UTRs","promoters","1to5kb","intergenic")))
g=ggplot(tmp4,aes(x=category,y=value,fill=annot))+geom_bar(position="stack",stat="identity",colour="white")+coord_flip()+theme_bw()+ylab("Frequency (%)")+xlab("")+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=10))+theme_classic()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),strip.background=element_blank(),strip.text.x=element_text(size=15))+facet_wrap(~clock)
ggsave("GenicAnnotationsByContribution+legend_25percent.pdf",g,width=8,height=8)
g=g+theme(legend.position="NONE")
ggsave("GenicAnnotationsByContribution_25percent.pdf",g,width=8,height=1.5)

# calc/get coefficient of variation and other stats
tmp1=read.table("path/to/epiAgesStats.tsv",header=T,sep="\t",stringsAsFactors=F)
tmp1$cv=tmp1$sd/tmp1$mean
tmp1$clock2=sub("Levine2018","PhenoAge clock",sub("Horvath2018","Skin & blood clock",sub("Horvath2013","Pan-tissue clock",tmp1$clock)))
write.table(tmp1,"epiAgesStats_mod.tsv",row.names=F,quote=F,sep="\t")








