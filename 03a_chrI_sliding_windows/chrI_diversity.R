## CHROMOSOME 1 NUCLEOTIDE DIVERSITY

library(ggplot2)
library("coin")

getwd()

chrI <- read.table("chrI.div",sep="\t")

head(chrI)

chrI$pi<-as.numeric(chrI$pi)
chrI$tajd<-as.numeric(chrI$tajd)

colnames(chrI)<-c("chr","window","n_SNPs","cov","pi","pop","tajd")

chrI$Population<-ifelse(chrI$pop=="SRR7470095","Marine","Freshwater")
chrI$Population <- factor(chrI$Population,levels=c("Marine","Freshwater"))

# plot code
piplotpop<-ggplot(chrI, aes(x = window, y = pi, colour = Population, group="Population")) +
  theme_bw()+
  labs(x="position on chrI",y="Pi of 10kb window",title = "A")+
  geom_line()

tajplotpop<-ggplot(chrI, aes(x = window, y = tajd, colour = Population, group="Population")) +
  theme_bw()+
  labs(x="position on chrI",y="Tajima's D of 10kb window",title = "A")+
  geom_line()

nrow(chrI)/2

piplotpop2<-ggplot(chrI, aes(x = Population, y = pi, fill = Population)) +
  theme_bw()+
  geom_boxplot() +
  labs(x="Population",y="Pi of 10kb window",title = "A")

tajplotpop2<-ggplot(chrI, aes(x = Population, y = tajd, fill = Population)) +
  theme_bw()+
  geom_boxplot() +
  labs(x="Population",y="Tajima's D of 10kb window",title = "B")

independence_test(pi~Population,data=chrI,distribution="approximate",alternative="greater")

independence_test(tajd~Population,data=chrI,distribution="approximate",alternative="less")
