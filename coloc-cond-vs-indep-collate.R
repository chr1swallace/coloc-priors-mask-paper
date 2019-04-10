#!/usr/bin/env Rscript
library(randomFunctions)
library(snpStats)
library(magrittr)
library(data.table)
devtools::load_all("~/RP/coloc")
source("~/DIRS.txt")

d <- COLOCINDEP

files <- list.files(d,full=TRUE,pattern="^csim.*.RData")
message("files found: ",length(files))
RESULTS <- vector("list",length(files))
for(i in seq_along(files)) {
    (load(files[i]))
    results[,nmeth:=.N,by=c("isim","method")]
    results[,nok:=all(nmeth==2),by="isim"]
    results[,f:=basename(files[i])]
    RESULTS[[i]] <- results
}
results <- rbindlist(RESULTS)
setnames(results,
         grep("PP",names(results),value=TRUE),
         gsub("PP.|.abf","",grep("PP",names(results),value=TRUE)))

## used r instead of r2
results[,r2.cvA1.hits:=r2.cvA1.hits^2]
results[,r2.cvB1.hits:=r2.cvB1.hits^2]
results[,r2.cvA2.hits:=r2.cvA2.hits^2]
results[,r2.cvB2.hits:=r2.cvB2.hits^2]

## try and label different signals according to max r2
results[,tested.cv1:="?"][,tested.cv2:="?"]
results[
r2.cvA1.hits > r2.cvB1.hits &r2.cvA1.hits >0.5,tested.cv1:="A"][
r2.cvB1.hits > r2.cvA1.hits &r2.cvB1.hits >0.5,tested.cv1:="B"][
r2.cvA2.hits > r2.cvB2.hits &r2.cvA2.hits >0.5,tested.cv2:="C"][
r2.cvB2.hits > r2.cvA2.hits &r2.cvB2.hits >0.5,tested.cv2:="D"]

results[,tested.cv1:=ifelse(r2.cvA1.hits > r2.cvB1.hits, "A", "B")]
results[,tested.cv2:=ifelse(r2.cvA2.hits > r2.cvB2.hits, "C", "D")]

## results[,tested.cv1:=0]
## results[,tested.cv1:=1]
## results[r2.tr1.cv1+r2.tr1.cv2>0.5,tested.cv1:=ifelse(r2.tr1.cv1>r2.tr1.cv2,1,2)]
## results[r2.tr2.cv1+r2.tr2.cv2>0.5,tested.cv2:=ifelse(r2.tr2.cv1>r2.tr2.cv2,1,2)]
results[,tested.cv:=paste(tested.cv1,tested.cv2,sep=":")]

## without multiple cv testing, just take first stepwise
best <- results[results[method=="step",.I[1],by=c("f","isim")]$V1,]
best[,method:="single"]
best[,tested.cv:="--"]
results <- rbind(results,best)

## simulations options
## 5 = share 1, weakest effect for one, strongest for other
## 4 = share 2, opposite effects
## 3 = share 1, equal effect (weakest) + indep each
## 2 = share 2, equal effects
## 1 = share 1, equal effect (strongest), + indep each
## indep want H4=1, H3=1, H1=1, H2=1
## step want H4=1, H3=3?
## 0 = share 0

sresults <- results[,.(h0=sum(H0),
                       h1=sum(H1) + sum(H2),
                       h3=sum(H3),
                       h4=sum(H4)),by=c("f","isim","method","special")]
m <- melt(sresults[,.(special,method,h0,h1,h3,h4)],
          c("special","method"))
m[,pc:=value/sum(value),by=c("special","method")]
  library(viridis)
library(ggplot2)
library(ggpubr)
ggplot(m,aes(x=method,y=pc,fill=variable)) + geom_col() +
  facet_grid(special ~.) + theme_pubr() +
scale_fill_viridis_d()
## next: look at best1, best2 and the CVs in details.
## 1. are the bests in LD with the CVs?
## 2. are

kk <- results[r2.tr1<0.01 & r2.tr2<0.01, .(h0=mean(H0),h1=mean(H1),h2=mean(H2),h3=mean(H3),h4=mean(H4)), by=c("special","method","tested.cv")]

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

ggsave("~/fig-twotwo.pdf",height=8,width=6)





results[, .(h0=sum(H0>0.5),h1=sum(H1>0.5),h2=sum(H2>0.5),h3=sum(H3>0.5),h4=sum(H4>0.5)), by=c("special","method")]
table(results$special)

## cv1 should be shared
summ <- results[, .(h0=sum(H0>0.5),h1=sum(H1>0.5),h2=sum(H2>0.5),h3=sum(H3>0.5),h4=sum(H4>0.5)), by=c("special","method","tested.cv")]
m <- melt(summ[,.(special,method,tested.cv,h0,h1,h2,h3,h4)],
          c("special","method","tested.cv"))
m[,pc:=value/sum(value),by=c("special","method","tested.cv")]

## distribution of best guess
cols <- c("h0","h1","h2","h3","h4")
results[,irow:=1:.N,by=c("f","isim","special","method")]
m <- melt(results,c("f","isim","special","method","irow"),toupper(cols))
m <- m[order(value,decreasing=TRUE)]
m <- m[ m[,.I[1],by=c("f","isim","method","irow")]$V1 ]
m[,value:=1]
m[,N:=.N,by=c("method","special")]
m <- m[,.(pc=value/N),by=c("variable","method","special")]
m[method!="single",pc:=4*pc]


library(viridis)
library(ggplot2)
library(ggpubr)
ggplot(m,aes(x=variable,y=pc,fill=variable)) + geom_col() +
  facet_grid(special ~ method,space="free_x",scales="free_x") + theme_pubr() +
scale_fill_viridis_d() + ggtitle("Average posterior in each comparison")
## next: look at best1, best2 and the CVs in details.
## 1. are the bests in LD with the CVs?
## 2. are

## average posterior
summ <- results[, .(h0=mean(H0),h1=mean(H1),h2=mean(H2),h3=mean(H3),h4=mean(H4)), by=c("special","method","tested.cv")]
m <- melt(summ[,.(special,method,tested.cv,h0,h1,h2,h3,h4)],
          c("special","method","tested.cv"))
m[,pc:=value,by=c("special","method","tested.cv")]


library(viridis)
library(ggplot2)
library(ggpubr)
ggplot(m,aes(x=tested.cv,y=pc,fill=variable)) + geom_col() +
  facet_grid(special ~ method,space="free_x",scales="free_x") + theme_pubr() +
scale_fill_viridis_d() + ggtitle("Average posterior in each comparison")
## next: look at best1, best2 and the CVs in details.
## 1. are the bests in LD with the CVs?
## 2. are

## what is distribution of expected H4 for each system?
cols <- c("h0","h1","h2","h3","h4")
summ <- copy(results)
setnames(summ,toupper(cols),cols)
summ <- summ[ , (cols) := lapply(.SD, sum), .SDcols = cols,by=c("f","isim","method")]
summ <- unique(summ,by=c("f","isim","method"))
m <- melt(summ, c("method","special"),cols)
library(viridis)
library(ggplot2)
library(ggpubr)
library(cowplot)
ggplot(m[method!="step",],aes(x=variable,y=value,fill=variable)) +
  ## geom_violin() +
  geom_point() + 
  geom_boxplot() +
  facet_grid(special ~ method,space="free_x",scales="free_x") + theme_pubr() +
scale_fill_viridis_d() + background_grid() + ggtitle("Sum of posterior probs in each analysis")

## what is expected H4 for each system?
summ <- results[, .(h4=mean(H4)), by=c("special","method")]
summ[method!="single",h4:=4*h4]
m <- melt(summ[,.(special,method,h4)],
          c("special","method"))
m[,pc:=value,by=c("special","method")]


library(viridis)
library(ggplot2)
library(ggpubr)
library(cowplot)
ggplot(m,aes(x=method,y=pc,fill=variable)) + geom_col() +
  facet_grid(special ~ .,space="free_x",scales="free_x") + theme_pubr() +
scale_fill_viridis_d() + background_grid()


library(ggplot2)
library(ggpubr)
ggplot(results,aes(x=cpp,y=ipp,col=rCV)) + geom_point() + geom_abline() + theme_pubr() + facet_wrap(ocv~beta)

results[ipp <0.5 & cpp>0.5,]

