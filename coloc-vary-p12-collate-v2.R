#!/usr/bin/env Rscript
library(randomFunctions)
library(snpStats)
library(magrittr)
library(data.table)
devtools::load_all("~/RP/coloc")
source("~/DIRS.txt")
args <- getArgs(defaults=list(N=1000,NSIM=100,NCV=3,SPECIAL=0),
                numeric=c("N","NSIM","NCV"))

d <- COLOCINDEP

byvars <- c("special","p12","N")
patt <- "cvaryp12-v4"
files <- list.files(d,pattern=patt,full=TRUE)
message("files found: ",length(files))
data <- lapply(files, function(f) {
    x=eval(as.symbol(load(f)))
    ## x[,f:=basename(f)]
    x})  %>% rbindlist(.,fill=TRUE)
data <- data[p12==1e-6,] ## adjust later
data[,PP.H4.abf:=unlist(PP.H4.abf)]
data[,p12:=as.numeric(unlist(p12))]
data[,nsnps:=as.numeric(unlist(nsnps))]
data[,PP.H0.abf:=unlist(PP.H0.abf)]
data[,PP.H1.abf:=unlist(PP.H1.abf)]
data[,PP.H2.abf:=unlist(PP.H2.abf)]
data[,PP.H3.abf:=unlist(PP.H3.abf)]
data[,PP.H4.abf:=unlist(PP.H4.abf)]
data[,Nsims:=.N,by=byvars]

## special=0..4 equiv to H_i

## probs <- function(nsnps,p12,p1=1e-4,p2=1e-4) {
##     h1 <- p1 * nsnps
##     h2 <- p2 * nsnps
##     h3 <- p1 * p2 * nsnps * (nsnps-1)
##     h4 <- p12 * nsnps
##     h0 <- 1 - h1 - h2 - h3 - h4
##     list(h0=h0,h1=h1,h2=h2,h3=h3,h4=h4)
## }
## data[,c("h0","h1","h2","h3","h4"):=probs(nsnps,p12)]
f <- function(p12)
    prior.snp2hyp(data$nsnps,p12=p12,p1=1e-4,p2=1e-4)
summary(data$p12)
pr0 <- f(1e-6) # value used
D <- data[,.(PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf)]/pr0

newdata <- lapply(c(1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8), function(pp) {
    newpp <- D * f(pp)
    cbind(newpp/rowSums(newpp),p12=pp,data[,.(N,special,ld,pop,NSNP)])
})  %>% rbindlist()

head(newdata)
tail(newdata)

## Jamie's question
library(ggplot2)
tmp=newdata[special==3,]
ggplot(tmp, aes(x=factor(p12),y=PP.H4.abf)) + geom_boxplot()

tmp[,.(prob.gt.0.9=mean(PP.H4.abf>0.9)),by=c("p12","N")]

library(viridis)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(scales)


cols <- grep("PP",names(newdata),value=TRUE)
ypp <- newdata[ , .(H0=mean(PP.H0.abf),
                   H1=mean(PP.H1.abf),
                   H2=mean(PP.H2.abf),
                   H3=mean(PP.H3.abf),
                   H4=mean(PP.H4.abf)), by=byvars]
y <- melt(ypp,byvars,value.name = "post",variable.name="hyp")
y[,x:=as.numeric(factor(N))]
tt <- unique(y[,.(x,N)])
tt[order(x),]
setnames(y,"special","sim")
y[,sim:=paste0("H",sim)]

pr <- sort(unique(newdata$p12))
pr
lab12 <- c("1 %*% 10^{-8}","5 %*% 10^{-8}",
              "1 %*% 10^{-7}","5 %*% 10^{-7}",
              "1 %*% 10^{-6}","5 %*% 10^{-6}",
              "1 %*% 10^{-5}","5 %*% 10^{-5}",
           "1 %*% 10^{-4}")  %>% paste("p[12] ==",.)

 
y[,w:=p12]
y[,p12:=NULL]
for(i in seq_along(pr))
    y[w==pr[[i]],p12:=lab12[[i]]]
y[,sim:=paste0("sim: ",sim)]
table(y$p12)
y[,p12:=factor(p12,levels=lab12)]

ggplot(y,aes(x=x,y=post,fill=factor(hyp),col=factor(hyp))) +
  geom_col(position="stack") +
  facet_grid(sim~p12,labeller=label_parsed) +
  scale_x_continuous("Sample size",breaks=c(1,4,7),
                     labels=c(expression(10^{2}), expression(10^{3}),expression(10^{4}))) +
  scale_y_continuous("Posterior probability",breaks=c(0,0.5,1),labels=c("","0.5","1")) +
  scale_fill_viridis_d("Posterior hypothesis") +
  scale_colour_viridis_d("Posterior hypothesis") +
  theme_pubr() +
  theme(axis.text=element_text(size=8),
        strip.text.x=element_text(size=8,face="bold"),
        strip.text.y=element_text(size=10,face="bold"),
        strip.background = element_rect(linetype="blank",fill="white"))
  ## scale_x_log10()

ggsave("~/fig-coloc-vary-p12.pdf",height=6*8/9,width=8)
