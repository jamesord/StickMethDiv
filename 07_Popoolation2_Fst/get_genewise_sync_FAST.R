args = commandArgs(trailingOnly=TRUE)

sync<-read.table(args[1])
sites<-read.table(args[2])

merged<-merge(sites,sync,by=c("V1","V2"))

newdata<-NULL

for (i in levels(as.factor(merged$V3.x))){
  cat<-subset(merged,V3.x==i)
  cat$rownum <- seq.int(nrow(cat))
  newdata<-rbind(newdata,cat)
  }

write.table(newdata[c(3,7,4,5,6)],file=args[3],quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")