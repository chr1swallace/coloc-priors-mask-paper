#!/usr/bin/env Rscript
library(randomFunctions)
library(snpStats)
library(magrittr)
library(data.table)
devtools::load_all("~/RP/coloc")
source("~/DIRS.txt")

d <- COLOCINDEP

files <- list.files(d,full=TRUE,pattern="^c12sim.*.RData")
message("files found: ",length(files))
RESULTS <- vector("list",length(files))
for(i in seq_along(files)) {
    (load(files[i]))
    if(!nrow(results))
        next
    results[,isim:=i]
    results[,nmeth:=.N,by=c("isim","method")]
    ## results[,nok:=all(nmeth==2),by="isim"]
    results[,f:=basename(files[i])]
    RESULTS[[i]] <- results
}
results <- rbindlist(RESULTS)
table(results$special)
setnames(results,
         grep("PP",names(results),value=TRUE),
         gsub("PP.|.abf","",grep("PP",names(results),value=TRUE)))

## used r instead of r2
results[,r2.cvA.hits:=r2.cvA.hits^2]
results[,r2.cvB.hits:=r2.cvB.hits^2]

## simulations options
## 5 = share 1, weakest effect for one, strongest for other
## 4 = share 2, opposite effects
## 3 = share 1, equal effect (weakest) + indep each
## 2 = share 2, equal effects
## 1 = share 1, equal effect (strongest), + indep each
## indep want H4=1, H3=1, H1=1, H2=1
## step want H4=1, H3=3?
## 0 = share 0

## without multiple cv testing, just take first stepwise
results[,tested.cv:=ifelse(r2.cvA.hits > r2.cvB.hits, "A", "B")]
best <- results[results[method=="step",.I[1],by=c("f","isim")]$V1,]
best[,method:="single"]
best[,tested.cv:="--"]
results <- rbind(results,best)

results[,nok:=TRUE]
results[nok==TRUE, .(h0=sum(H0>0.5),h1=sum(H1>0.5),h2=sum(H2>0.5),h3=sum(H3>0.5),h4=sum(H4>0.5)), by=c("special","method")]

## results[,tested.cv1:=0]
## results[r2.tr1.cv1+r2.tr1.cv2>0.5,tested.cv1:=ifelse(r2.tr1.cv1>r2.tr1.cv2,1,2)]
## cv1 should be shared
## ss <- strsplit(results$details,"/")
## l <- sapply(ss,length)
## table(l)

## cv1 <- sapply(ss,"[[",1)
## cv2 <- sapply(ss,"[[",2)
## cv3 <- sapply(ss,"[[",5)

kk <- results[, .(h0=mean(H0),h1=mean(H1),h2=mean(H2),h3=mean(H3),h4=mean(H4)), by=c("special","method","tested.cv")]

m <- melt(kk,c("special","method","tested.cv"))
m[,method:=factor(method,levels=c("single","cond","mask"))]
m[,special:=factor(special,levels=c(0,1,3,5,2,4))]
## levels(m$special) <- c("none","one, strong","one,weak")
ggplot(m[method!="step",],aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special ~ method,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hypoth.") +
  ylab("Average posterior in each comparison") +
  xlab("Tested variants") +
  theme(legend.position="right")

ggsave("~/fig-onetwo.pdf",height=8*4/6,width=6)


summ <- results[nok==TRUE, .(h0=sum(H0>0.5),h1=sum(H1>0.5),h2=sum(H2>0.5),h3=sum(H3>0.5),h4=sum(H4>0.5)), by=c("special","method","tested.cv")]

summ[,h13:=h1+h3]
m <- melt(summ[,.(special,method,tested.cv,h0,h1,h2,h3,h4)],
          c("special","method","tested.cv"))
m[,pc:=value/sum(value),by=c("special","method","tested.cv")]

library(viridis)
library(ggplot2)
library(ggpubr)
ggplot(m[method!="step",],aes(x=tested.cv,y=pc,fill=variable)) + geom_col() +
  facet_grid(special ~ method,scales="free_x",space="free_x") + theme_pubr() +
scale_fill_viridis_d()
## next: look at best1, best2 and the CVs in details.
## 1. are the bests in LD with the CVs?
## 2. are


## results[,cor(cpp,ipp),by=c("ocv","beta")]

## library(ggplot2)
## library(ggpubr)
## ggplot(results,aes(x=cpp,y=ipp,col=rCV)) + geom_point() + geom_abline() + theme_pubr() + facet_wrap(ocv~beta)

## results[ipp <0.5 & cpp>0.5,]

