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

byvars <- c("special","p12","N","ld","pop")
patt <- "cvaryp12"
files <- list.files(d,pattern=patt,full=TRUE)
message("files found: ",length(files))
data <- lapply(files, function(f) {
    x=eval(as.symbol(load(f)))
    ## x[,f:=basename(f)]
    x})  %>% rbindlist(.,fill=TRUE)
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

weights <- function(special,nsnps,p12,p1=1e-4,p2=1e-4) {
    h1 <- p1 * nsnps
    h2 <- p2 * nsnps
    h3 <- p1 * p2 * nsnps * (nsnps-1)
    h4 <- p12 * nsnps
    h0 <- 1 - h1 - h2 - h3 - h4
    h <- cbind(h0,h1,h2,h3,h4)
    h[ cbind(1:nrow(h) , as.numeric(special)+1) ] 
}
for(ip in unique(data$p12)) {
    thing <- paste0("w.",ip)
    data[,(thing):=weights(special,nsnps,ip)]
}
## data[,w:=NULL]
## setnames(data,"w.5e-06","w")
## data[,w:=weights(N,special,nsnps,p12)]

## data[,right:=ifelse(special=="0",PP.H0.abf,
##              ifelse(special=="1",PP.H1.abf,
##              ifelse(special=="2",PP.H2.abf,
##              ifelse(special=="3",PP.H3.abf,
##                              PP.H4.abf))))]
## data[,right0:=right/(1-PP.H0.abf)]

library(ggplot2)
x <- melt(data[,c("Nsims",byvars,
                  grep("PP.",names(data),value=TRUE),
                  grep("^w",names(data),value=TRUE),"p12"),with=FALSE],
          id.vars=c("Nsims",byvars,
                    grep("PP.",names(data),value=TRUE)),
          measure.vars=grep("^w",names(data),value=TRUE),
          variable.name="weight.p12",
          value.name="weight")


## prior
ypr <- x[,list(pr=sum(weight/Nsims)),by=byvars]
## tt <- table(ypr$nsnps)
sort(tt,decreasing=TRUE)
ypr <- unique(ypr[,.(special,p12,pr)])
ypr[,sum(pr),by=c("p12")] # should be all 1

library(viridis)
ypr[,x:=1]
ggplot(ypr,aes(x=x,y=pr,fill=special)) + geom_col() + scale_fill_viridis_d("post. hyp")  + facet_grid(.~p12) + ggtitle("per SNP priors")


## posterior
x[,special:=NULL]
x <- melt(x,
          id.vars=c("Nsims","weight",setdiff(byvars,"special")),
          measure.vars = 
          variable.name="hyp",
          value.name="pp")
## x[,special:=gsub("PP.H|.abf","",hyp)]
x[,special:=substr(hyp,5,5)] # always character 5, substr is faster
ypp <- x[,list(pp=sum(pp*weight/Nsims)),by=c("p12",byvars)]
ypp[,sum(pp),by=c("weight.p12","p12",setdiff(byvars,"special"))] # should be all 1



y <- merge(ypr,ypp,by=c("weight.p12","special"))
tmp=y[weight.p12=="w.1e-04" & p12 %in% c(1e-04,1e-08),]
tmp

## ggplot(y[weight.p12=="w.1e-04",],aes(x=pr,y=pp)) + geom_abline() + geom_point(aes(colour=special)) + facet_wrap(~p12)

## ggplot(y,aes(x=pr,y=pp)) + geom_abline() + geom_point(aes(colour=special)) + facet_grid(weight.p12~p12)

y <- x[special!="0",list(mv=sum(right*value)/sum(value)),by=c("variable","p12")]
ggplot(y,aes(x=variable,y=mv,col=factor(p12),group=factor(p12))) + geom_point() + geom_path() + facet_wrap(~p12)


y0 <- x[special!="0",list(mv=sum(right0*value)/sum(value)),by=c("variable","p12")]
ggplot(y0,aes(x=p12,y=mv,col=variable)) + geom_point() + geom_path()




y0 <- data[special!="0",list(mv=sum(right0*w)/sum(w)),
          by=c("p12")]

ggplot(y,aes(x=p12,y=mv)) + geom_point() + ylim(0.5,1)
ggplot(y0,aes(x=p12,y=mv)) + geom_point() + ylim(0.5,1)

y <- data[,.(mv=mean(right)),by=c("p12")]
ggplot(y,aes(x=p12,y=mv)) + geom_point() + ylim(0.5,1)

y <- x[,
x <- melt(data[,c("right","special","p12","w"),with=FALSE],
          c("special","p12","w"))
x[,variable:=gsub("PP.|.abf","",variable)]
x[,special:=as.integer(special)]
x[,p12:=as.factor(unlist(p12))]
x[,value:=unlist(value)]


## y <- x[,.(mv=mean(value)),by=c("special","p12","variable")]


library(ggplot2)
x <- melt(data[,c(grep("PP",names(data),value=TRUE),"special","p12","w"),with=FALSE],
          c("special","p12","w"))
x[,variable:=gsub("PP.|.abf","",variable)]
x[,special:=as.integer(special)]
x[,p12:=as.factor(unlist(p12))]
x[,value:=unlist(value)]


## y <- x[,.(mv=mean(value)),by=c("special","p12","variable")]
y <- x[,.(mv=sum(value*w)/sum(w)),by=c("special","p12","variable")]
ggplot(y,aes(x=special,y=mv,fill=variable)) + geom_col() + facet_wrap(~p12)



data[,
y <- data[,.(mv=mean(PP.H4.abf)),by=c("p12")]
