# characters related to sample ID are masked by "__"

#-------------------
# load packages
#-------------------
print("load packages")
library(reshape2)
library(ggplot2)

#-------------------
# load datasets
#-------------------
serum=read.table("path/to/serum/test/results",head=T,stringsAsFactors=F)
serum1=serum;serum1$sample=sub("__","__",sub("__","__",serum1$sample))
serum2=serum;serum2$sample=sub("__","__",sub("__","__",serum2$sample))
serumMod=rbind(serum1,serum2)
epiage=read.table("path/to/epiAges.tsv",head=T,stringsAsFactors=F,sep="\t")
cellcount1=read.table("path/to/cellcout/personA",head=T,stringsAsFactors=F)
cellcount2=read.table("path/to/cellcout/personB",head=T,stringsAsFactors=F)
cellcount=rbind(cellcount1,cellcount2)
cellcount$id=rownames(cellcount)

#-------------------
# plot
#-------------------
tmp1=cellcount
tmp1$person="Person A"
tmp1[grep("__|__",tmp1$id),"person"]="Person B"
tmp1$cell="PBMC"
tmp1[grep("__|__",tmp1$id),"cell"]="Monocyte"
tmp1$day=as.numeric(sub("^[A-Z][A-Z]","",tmp1$id))
tmp2=tmp1[,c("person","day","cell","CD8T","CD4T","NK","Bcell","Mono","Gran")]
tmp3=melt(tmp2,id.vars=c("person","day","cell"))
g=ggplot(tmp3,aes(x=day,y=value,colour=variable))+geom_line()+geom_point(shape=16,size=0.3)+theme_bw()+facet_grid(cell~person)+xlab("Day")+ylab("Proportion")+scale_y_continuous(breaks=seq(0,1,0.2),limits=c(-0.05,1.05))
ggsave("cellCount+legend.pdf",g,width=4,height=4)
g=g+theme(legend.position="NONE")
ggsave("cellCount.pdf",g,height=4,width=4)


#-------------------
# age vs cellcount (consider only PBMC)
#-------------------
tmp1=NULL
tmp3=NULL
for(clock in unique(epiage$clock)){
 for(type in unique(epiage$type)){
  tmpA=epiage[epiage$clock==clock & epiage$type==type & epiage$cell=="PBMC",]
  tmpB=merge(tmpA[,c("id","EpiAge","person")],cellcount[,c("id","CD8T","CD4T","NK","Bcell","Mono","Gran")],by="id",all.x=T,all.y=F)
  tmpC=lm("EpiAge~CD8T+CD4T+NK+Bcell+Mono",data=tmpB)
  tmpB$adjEpiAge=as.numeric(summary(tmpC)$residuals)
  tmpB$day=as.numeric(sub("^[A-Z][A-Z]","",tmpB$id))
  tmpC=tmpB[,c("id","person","EpiAge","adjEpiAge","day")]
  tmpC$clock=clock
  tmpC$type=type
  tmp3=rbind(tmp3,tmpC)
  tmpB1=tmpB[tmpB$person=="Person A",]
  tmpB2=tmpB[tmpB$person=="Person B",]
  sdRaw1=sd(tmpB1$EpiAge)
  sdAdj1=sd(tmpB1$adjEpiAge)
  ftest1=var.test(tmpB1$EpiAge,tmpB1$adjEpiAge)$p.value
  sdRaw2=sd(tmpB2$EpiAge)
  sdAdj2=sd(tmpB2$adjEpiAge)
  ftest2=var.test(tmpB2$EpiAge,tmpB2$adjEpiAge)$p.value
  range1=max(tmpB1$adjEpiAge)-min(tmpB1$adjEpiAge)
  range2=max(tmpB2$adjEpiAge)-min(tmpB2$adjEpiAge)
  tmpD=NULL
  for(n in 2:24){
   tmpB1=tmpB1[order(tmpB1$day,decreasing=F),]
   tmpB2=tmpB2[order(tmpB2$day,decreasing=F),]
   tmpC1=tmpB1[c((n-1):n),]
   tmpC2=tmpB1[c((n-1):n),]
   tmpD1=abs(as.numeric(tmpC1$adjEpiAge)[1]-as.numeric(tmpC1$adjEpiAge)[2])/(as.numeric(tmpC1$day)[2]-as.numeric(tmpC1$day)[1])
   tmpD2=abs(as.numeric(tmpC2$adjEpiAge)[1]-as.numeric(tmpC2$adjEpiAge)[2])/(as.numeric(tmpC2$day)[2]-as.numeric(tmpC2$day)[1])
   day1=tmpC1[1,"day"]
   day2=tmpC1[2,"day"]
   tmpD=rbind(tmpD,c(day1,day2,tmpD1,tmpD2))
  }
  daily1=max(tmpD[,3])
  daily2=max(tmpD[,4])
  res=c(clock,type,range1,sdRaw1,sdAdj1,ftest1,daily1,range2,sdRaw2,sdAdj2,ftest2,daily2)
  tmp1=rbind(tmp1,res)
 }
}
tmp2=as.data.frame(tmp1,stringsAsFactors=F)
colnames(tmp2)=c("clock","type","range1","sdRaw1","sdAdj1","ftest1","daily1","range2","sdRaw2","sdAdj2","ftest2","daily2")
write.table(tmp2,"adjustedEpiAgesSummary.tsv",row.names=F,quote=F,sep="\t")
ages=tmp3

#-------------------
# age vs serum test results (RawEpiage of Mono/PBMC and AdjEpiAge of PBMC)
#-------------------
ages$cell="PBMC"
tmp1=epiage[epiage$cell=="Monocyte",c("id","person","EpiAge","day","clock","type","cell")]
tmp1$adjEpiAge=NA
tmp2=rbind(ages,tmp1)
tmp3=NULL
for(clock in unique(epiage$clock)){
 for(type in unique(epiage$type)){
  for(data in c("mono","pbmc","adjPbmc")){
   if(data=="mono"){
    tmpA=tmp2[tmp2$clock==clock & tmp2$type==type & tmp2$cell=="Monocyte",c("id","person","EpiAge")]
   }else if(data=="pbmc"){
    tmpA=tmp2[tmp2$clock==clock & tmp2$type==type & tmp2$cell=="PBMC",c("id","person","EpiAge")]
   }else if(data=="adjPbmc"){
    tmpA=tmp2[tmp2$clock==clock & tmp2$type==type & tmp2$cell=="PBMC",c("id","person","adjEpiAge")]
    colnames(tmpA)[3]="EpiAge"
   }
   tmpB=merge(tmpA,serumMod,by.x="id",by.y="sample",all=F)
   for(variable in colnames(serumMod)[-3:-1]){
    tmpC=tmpB[,c("person","EpiAge",variable)]
    colnames(tmpC)[3]="variable"
    lm1=lm(EpiAge~variable,data=tmpC[tmpC$person=="Person A",])
    lm2=lm(EpiAge~variable,data=tmpC[tmpC$person=="Person B",])
    coef1=as.numeric(summary(lm1)$coefficients[2,1])
    se1=as.numeric(summary(lm1)$coefficients[2,2])
    pval1=as.numeric(summary(lm1)$coefficients[2,4])
    coef2=as.numeric(summary(lm2)$coefficients[2,1])
    se2=as.numeric(summary(lm2)$coefficients[2,2])
    pval2=as.numeric(summary(lm2)$coefficients[2,4])
    res=c(clock,type,data,variable,coef1,se1,pval1,coef2,se2,pval2)
    tmp3=rbind(tmp3,res)
   }
  }
 }
}
tmp4=as.data.frame(tmp3,stringsAsFactors=F)
colnames(tmp4)=c("clock","type","data","variable","coef1","se1","pval1","coef2","se2","pval2")
# insert R codes for meta-analysis (provided by Dr. Matti Pirinen (https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html, accessed on 6 April 2022)) and make "tmp5" including meta-analysis results
write.table(tmp5,"EpiageSerumAssoc.tsv",row.names=F,quote=F,sep="\t")












