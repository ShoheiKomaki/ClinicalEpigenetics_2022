# characters related to sample ID are masked by "__"

# load packages
library(ggplot2)
library(data.table)
library(reshape2)

# load datasets (lists of clock CpGs and their coefficient + intercept)
horvath2013=read.csv("path/to/file",header=T,stringsAsFactors=F)
horvath2018=read.csv("path/to/file",header=T,stringsAsFactors=F)
levine2018=read.csv("path/to/file",header=T,stringsAsFactors=F)

tmp00=NULL
for(type in c("raw","bmiq","horvath")){
 print(paste0("normalization type: ",type))
 print("loading dataset")
 if(type=="raw"){
  type2="Raw data"
  fn="path/to/raw/data"
  tmp1=as.data.frame(fread(fn),stringsAsFactors=F)
  colnames(tmp1)[1]="cgID"
 }else if(type=="bmiq"){
  type2="BMIQ normalization"
  fn="path/to/bmiq/data"
  tmp1=as.data.frame(fread(fn),stringsAsFactors=F)
  colnames(tmp1)[1]="cgID"
 }else if(type=="horvath"){
  type2="Horvath normalization"
  fn="path/to/horvath/data"
  tmp1=as.data.frame(fread(fn),stringsAsFactors=F)
  colnames(tmp1)[1]="cgID"
 }
 df1=tmp1
 print("performing PCA")
 for(target in c("allCpG","Horvath2013","Horvath2018","Levine2018")){
  if(target=="allCpG"){
   tmp1=df1
  }else if(target=="Horvath2013"){
   tmp1=df1[df1$cgID%in%horvath2013$CpGmarker[-1],]
  }else if(target=="Horvath2018"){
   tmp1=df1[df1$cgID%in%horvath2018$ID[-1],]
  }else if(target=="Levine2018"){
   tmp1=df1[df1$cgID%in%levine2018$CpG[-1],]
  }
  tmp2=apply(tmp1,1,function(x){sd(x[-1],na.rm=T)})
  tmp3=tmp1[tmp2!=0,]
  pcares=prcomp(t(na.omit(tmp3[,-1])))
  p2=as.data.frame(pcares$x,stringsAsFactors=F)
  expPC1=summary(pcares)$importance[2,"PC1"]
  expPC2=summary(pcares)$importance[2,"PC2"]
  expPC3=summary(pcares)$importance[2,"PC3"]
  expPC4=summary(pcares)$importance[2,"PC4"]
  p3=p2
  p3$sample=NA
  p3[grep("__",rownames(p3)),"sample"]="PBMC (person B)"
  p3[grep("__",rownames(p3)),"sample"]="Monocyte (person B)"
  p3[grep("__",rownames(p3)),"sample"]="PBMC (person A)"
  p3[grep("__",rownames(p3)),"sample"]="Monocyte (person A)"
  g=ggplot()+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC1,y=PC2),colour="#159fd2")+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC1,y=PC2),colour="#ee6129")+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC1,y=PC2),colour="#2f5597")+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC1,y=PC2),colour="#c00000")+
   theme_bw()+
   xlab(paste0("PC1 (",100*expPC1,"%)"))+
   ylab(paste0("PC2 (",100*expPC2,"%)"))+
   theme(legend.position="NONE")
  ggsave(paste0(type,"_dataset_",target,"PC1-2.pdf"),g,height=2,width=2)
  g=ggplot()+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC3,y=PC4),colour="#159fd2")+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC3,y=PC4),colour="#ee6129")+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC3,y=PC4),colour="#2f5597")+
   geom_point(p3[grep("__",rownames(p3)),],
    mapping=aes(x=PC3,y=PC4),colour="#c00000")+
   theme_bw()+
   xlab(paste0("PC3 (",100*expPC3,"%)"))+
   ylab(paste0("PC4 (",100*expPC4,"%)"))+
   theme(legend.position="NONE")
  ggsave(paste0(type,"_dataset_",target,"PC3-4.pdf"),g,height=2,width=2)
 }
 print("calcurating epigenetic ages")
 tmp0=NULL
 day=c(1,3,8,10,14,17,22,24,29,31,32,38,42,43,46,53,56,57,72,77,78,79,NA,81,84)
 for(target in c("Horvath2013","Horvath2018","Levine2018")){
  if(target=="Horvath2013"){
   int=horvath2013[1,2]
   tmp1=horvath2013[-1,c("CpGmarker","CoefficientTraining")]
   colnames(tmp1)=c("cgID","coef")
   tmp2=merge(tmp1,df1,by="cgID")
  }else if(target=="Horvath2018"){
   int=horvath2018[1,2]
   tmp1=horvath2018[-1,c("ID","Coef")]
   colnames(tmp1)=c("cgID","coef")
   tmp2=merge(tmp1,df1,by="cgID")
  }else if(target=="Levine2018"){
   int=levine2018[1,6]
   tmp1=levine2018[-1,c("CpG","Weight")]
   colnames(tmp1)=c("cgID","coef")
   tmp2=merge(tmp1,df1,by="cgID")
  }
  tmp3=NULL
  for(id in colnames(tmp2)[-2:-1]){
   tmpA=tmp2[,c("coef",id)]
   tmpA$tmp=tmpA$coef * tmpA[,id]
   if(target!="Levine2018"){
    tmpB=21 * (sum(tmpA$tmp) + int) + 20
   }else{
    tmpB=sum(tmpA$tmp)+int
   }
   tmpC=nrow(tmpA[is.na(tmpA$tmp),])
   tmpD=day[(as.numeric(sub("^[A-Z][A-Z]","",id)))-11]
   tmp3=rbind(tmp3,c(id,tmpB,tmpC,tmpD))
  }
  tmp4=as.data.frame(tmp3,stringsAsFactors=F)
  colnames(tmp4)=c("id","EpiAge","Missing","day")
  tmp4$clock=target
  tmp0=rbind(tmp0,tmp4)
 }
 tmp0$type=type2
 tmp00=rbind(tmp00,tmp0)
}
tmp00$person=sub("__","Person B",sub("__","Person A",sub(".$","",sub("[0-9][0-9]$","",tmp00$id))))
tmp00$cell=sub("P","PBMC",sub("M","Monocyte",sub("^.","",sub("[0-9][0-9]$","",tmp00$id))))
df2=tmp00
write.table(df2,"epiAges.tsv",row.names=F,quote=F,sep="\t")

g=ggplot(df2,aes(x=as.numeric(day),y=as.numeric(EpiAge),colour=cell,linetype=type))+geom_line()+facet_grid(clock~person)+theme_bw()+ylab("DNA methylation age")+xlab("Day")
ggsave("epiAges+legend.pdf",g,height=4,width=4)
g=g+theme(legend.position="NONE")
ggsave("epiAges.pdf",g,height=4,width=4)

# get summary
df2=read.table("epiAges.tsv",header=T,sep="\t",stringsAsFactors=F)
tmp1=NULL
for(type in c("Raw data","BMIQ normalization","Horvath normalization")){
 for(clock in c("Horvath2013","Horvath2018","Levine2018")){
  for(cell in c("PBMC","Monocyte")){
   for(person in c("Person A","Person B")){
    tmpA=df2[df2$type==type & df2$clock==clock & df2$person==person & df2$cell==cell,]
    min=min(as.numeric(tmpA$EpiAge))
    max=max(as.numeric(tmpA$EpiAge))
    range=max-min
    mean=mean(as.numeric(tmpA$EpiAge))
    median=median(as.numeric(tmpA$EpiAge))
    sd=sd(as.numeric(tmpA$EpiAge))
    tmpD=NULL
    for(n in 2:nrow(tmpA)){
     tmpB=tmpA[c((n-1):n),]
     tmpC=abs(as.numeric(tmpB$EpiAge)[1]-as.numeric(tmpB$EpiAge)[2])/(as.numeric(tmpB$day)[2]-as.numeric(tmpB$day)[1])
     day1=tmpB[1,"day"]
     day2=tmpB[2,"day"]
     tmpD=rbind(tmpD,c(day1,day2,tmpC))
    }
    dayMax=tmpD[order(as.numeric(tmpD[,3]),decreasing=T)[1],]
    tmp1=rbind(tmp1,c(type,clock,cell,person,min,max,range,mean,median,sd,dayMax))
   }
  }
 }
}
tmp3=as.data.frame(tmp1,stringsAsFactors=F)
colnames(tmp3)=c("type","clock","cell","person","min","max","range","mean","median","sd","dayMaxFrom","dayMaxTo","dayMax")
df3=tmp3
write.table(df3,"epiAgesStats.tsv",row.names=F,quote=F,sep="\t")


# SD vs Coefficient relationships
tmp00=NULL
df4=NULL
for(type in c("raw","bmiq","horvath")){
 print(paste0("normalization type: ",type))
 print("loading dataset")
 if(type=="raw"){
  type2="Raw data"
  fn="path/to/raw/data"
  tmp1=as.data.frame(fread(fn),stringsAsFactors=F)
  colnames(tmp1)[1]="cgID"
 }else if(type=="bmiq"){
  type2="BMIQ normalization"
  fn="path/to/BMIQ/data"
  tmp1=as.data.frame(fread(fn),stringsAsFactors=F)
  colnames(tmp1)[1]="cgID"
 }else if(type=="horvath"){
  type2="Horvath normalization"
  fn="path/to/Horvath/data"
  tmp1=as.data.frame(fread(fn),stringsAsFactors=F)
  colnames(tmp1)[1]="cgID"
 }
 df1=tmp1
 print("preparing coefficients")
 tmp1=horvath2013[-1,c("CpGmarker","CoefficientTraining")]
 tmp1$clock="Horvath 2013"
 tmp2=horvath2018[-1,c("ID","Coef")]
 tmp2$clock="Horvath 2018"
 tmp3=levine2018[-1,c("CpG","Weight")]
 tmp3$clock="Levine 2018"
 colnames(tmp1)=colnames(tmp2)=colnames(tmp3)=c("cgID","coef","clock")
 tmp4=rbind(tmp1,tmp2,tmp3)
 tmp5=as.data.frame(t(apply(tmp4,1,function(x){
  cgID=as.character(x)[1]
  coef=as.numeric(x[2])
  clock=as.character(x[3])
  tmpA=df1[df1$cgID==cgID,]
  __=sd(as.numeric(tmpA[,grep("__",colnames(tmpA))]))
  __=sd(as.numeric(tmpA[,grep("__",colnames(tmpA))]))
  __=sd(as.numeric(tmpA[,grep("__",colnames(tmpA))]))
  __=sd(as.numeric(tmpA[,grep("__",colnames(tmpA))]))
  return(c(cgID,coef,clock,__,__,__,__))
 })),stringsAsFactors=F)
 colnames(tmp5)=c("cgID","coef","clock","AP","AM","BP","BM")
 tmp6=melt(tmp5,id.vars=c("cgID","coef","clock"))
 tmp6$cell="PBMC"
 tmp6[tmp6$variable=="AM" | tmp6$variable=="BM","cell"]="Monocyte"
 tmp6$person="Person A"
 tmp6[grep("B",tmp6$variable),"person"]="Person B"
 g=ggplot(tmp6,aes(y=as.numeric(coef),x=as.numeric(value),colour=cell))+
  geom_point(shape=16,size=0.5)+
  facet_grid(clock~person,scales="free")+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed")+
  ylab("Coefficient (weight)")+
  xlab("Standard deviation of DNA methylation level")
 ggsave(paste0("CpGEffect_",type,"+legend.pdf"),g,height=4,width=4)
 g=g+theme(legend.position="NONE")
 ggsave(paste0("CpGEffect_",type,".pdf"),g,height=4,width=4)
 tmp6$value2=as.numeric(tmp6$value)
 tmp6[tmp6$value2>0.06,"value2"]=0.06
 tmp6$shape="16"
 tmp6[tmp6$value2>0.06,"shape"]="18"
 g=ggplot(tmp6,aes(y=as.numeric(coef),x=as.numeric(value2),
   colour=cell,shape=shape))+
  geom_point(size=0.5)+
  facet_grid(clock~person,scales="free")+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed")+
  ylab("Coefficient (weight)")+
  xlab("Standard deviation of DNA methylation level")+
  theme(legend.position="NONE")+
  xlim(0,0.06)+
  scale_shape_manual(values=c("16"=16,"18"=18))
 ggsave(paste0("CpGEffect_",type,"_trim.pdf"),g,height=4,width=4)
 tmp7=apply(tmp5[,4:7],1,function(x){return(mean(as.numeric(x)))})
 tmp5$meanSD=as.numeric(tmp7)
 tmp5$type=type
 df4=rbind(df4,tmp5[,c("type","cgID","coef","clock","meanSD")])
}
write.table(df4,"SDvsCoef.tsv",row.names=F,quote=F,sep="\t")

# compare PBMC and monocyte ages
tmp1=read.table("epiAges.tsv",header=T,sep="\t",stringsAsFactors=F)
tmp2=NULL
for(type in unique(tmp1$type)){
 for(clock in unique(tmp1$clock)){
  for(person in unique(tmp1$person)){
   tmpA=tmp1[tmp1$type==type & tmp1$clock==clock & tmp1$person==person,]
   tmpPBMC=tmpA[tmpA$cell=="PBMC",c("day","EpiAge")]
   colnames(tmpPBMC)[2]="PBMC"
   tmpMono=tmpA[tmpA$cell=="Monocyte",c("day","EpiAge")]
   colnames(tmpMono)[2]="Mono"
   tmpB=merge(tmpPBMC,tmpMono,by="day")
   tmpC=t.test(tmpB$PBMC,tmpB$Mono,paired=T)
   tval=as.numeric(tmpC$statistic)
   pval=as.numeric(tmpC$p.value)
   tmpD=c(type,clock,person,tval,pval)
   tmp2=rbind(tmp2,tmpD)
  }
 }
}
tmp3=as.data.frame(tmp2,stringsAsFactors=F)
colnames(tmp3)=c("type","clock","person","t_value","p_value")
write.table(tmp3,"paired-t.tsv",row.names=F,quote=F,sep="\t")














