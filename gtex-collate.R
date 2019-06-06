#!/usr/bin/env Rscript

## d <- "/home/cew54/share/Data/gtex"
library(magrittr)
library(randomFunctions)
devtools::load_all("~/RP/GUESSFM")

args <- getArgs(defaults=list(chr=22))

## d <- "/home/cew54/share/Data/gtex/eur-Whole_Blood-001"
d <- "/home/cew54/scratch/gtex-gfm"


exp.ncv <- function(d,snps,dist) {
    snps.drop <- rownames(snps)[abs(snps$pdiff)>dist]
    snps.keep <- setdiff(rownames(snps),snps.drop)
    if(length(snps.keep) < length(snps.drop))
        modcount <- sapply(d@model.snps, function(s) sum(s %in% snps.keep))
    else
        modcount <- d@models$size - sapply(d@model.snps, function(s) sum(s %in% snps.drop))
    c(sum(modcount * d@models$PP),length(snps.keep))
}

if(!interactive()) {
    library(parallel)
    options(mc.cores=4)
    chr <- args$chr
    fdd <- paste0("~/scratch/dd-",chr,".RData")
    fsd <- paste0("~/scratch/snpdata-",chr,".RData")
    fgt <- paste0("~/scratch/gtex-",chr,".RData")
    ## if(file.exists(fgt))
    ##     next
    files <- list.files(file.path(d,paste0("chr",chr)),full=TRUE,pattern="^[0-9]+$")
    if(!file.exists(fdd)) {
        dd <- lapply(files,read.snpmod)
        save(dd,file=fdd)
    } else {
        (load(fdd))
    }

    if(!file.exists(fsd)) {
        snpdata <- lapply(files, function(dd) {
            (load(file.path(dd,"data.RData")))
            snps$pdiff <- snps$position - thisgene$start
            snps[,c("position","pdiff")]
        })
        save(snpdata,file=fsd)
    } else {
        (load(fsd))
    }

    en <- mclapply(c((1:20)*1e+5), function(dthr) {
        message(dthr)
        mapply(exp.ncv,dd,snpdata,dthr,SIMPLIFY=FALSE)
    }    )

    save(en,file=fgt)

    q("no")
}


## make plots
    fgt <- paste0("~/scratch/gtex-",1:22,".RData")
file.exists(fgt)


fgt <- fgt[ file.exists(fgt) ]
data <- lapply(fgt, function(f) {
    tmp <- eval(as.symbol(load(f))) %>% lapply(., function(x) {
        x <- do.call("rbind",x) 
    }) 
})
str(data)

## data <- data[[21]]
library(data.table)
library(ggplot2)
library(cowplot)
library(ggpubr)
ncv <- lapply(data,function(y) lapply(y, function(x) x[,1]))  %>% as.data.frame()
totalsnps <- lapply(data,function(y) lapply(y, function(x) x[,2]))  %>% as.data.frame()
fraccv <- lapply(data, function(y) lapply(y, function(x) x[,1]/x[,2] ))  %>% as.data.frame()
names(ncv) <- names(fraccv) <- sprintf("d-%012d",c((1:20)*1e+5))

f <- function(data) {
    data <- melt(data)
    data$dist <- as.numeric(sub("d-0*","",data$variable))
    as.data.table(data)
}
ncv %<>% f()
fraccv %<>% f()
data <- cbind(ncv,frac=fraccv$value)

y <- data[,.(qn=median(value,na.rm=TRUE),
             lqn=quantile(value,0.25,na.rm=TRUE),
             uqn=quantile(value,0.75,na.rm=TRUE),
             qf=median(frac,na.rm=TRUE),
             lqf=quantile(frac,0.25,na.rm=TRUE),
             uqf=quantile(frac,0.75,na.rm=TRUE) ),by="dist"]
y[,what:="GTeX whole blood eQTL"]


library(ggrepel)
p <- ggplot(y,aes(x=dist/1e+3,y=qn)) +
  geom_pointrange(aes(ymin=lqn,ymax=uqn)) +
  geom_path() + geom_hline(yintercept=1e-4,linetype="dashed",col="darkblue") +
  labs(x="Distance to TSS (kb)",y="Estimate of q") +
  ## geom_label_repel(aes(x=x,label=lab),
  ##                  data=data.frame(x=950,qn=1e-4,lab="q == 10^{-4}"),
  ##                  hjust=1.2,vjust=-3,#fontface="bold",
  ##                  col="darkblue",size=2.5,parse=TRUE) +
  facet_grid(. ~ what) +
  ## scale_y_log10() +
  scale_x_continuous(breaks=c(500,1000,1500,2000)) +
  theme_pubr() +
  background_grid() +
  theme(axis.text=element_text(size=8))
p

q <- ggplot(y,aes(x=dist/1e+3,y=qf)) +
  geom_pointrange(aes(ymin=lqf,ymax=uqf)) +
  geom_path() + geom_hline(yintercept=1e-4,linetype="dashed",col="darkblue") +
  labs(x="Distance to TSS (kb)",y="Estimate of q") +
  geom_label_repel(aes(x=x,label=lab),
                   data=data.frame(x=1950,qf=1e-4,lab="q == 10^{-4}"),
                   hjust=1.2,vjust=-3,#fontface="bold",
                   col="darkblue",size=2.5,parse=TRUE) +
  facet_grid(. ~ what) +
  ## scale_y_log10() +
  scale_x_continuous(breaks=c(500,1000,1500,2000)) +
  theme_pubr() +
  background_grid() +
  theme(axis.text=element_text(size=8))
## annotate(geom="text",x=50,y=1e-4,label="q == 10^{-4}",vjust=-1,col="darkblue",parse=TRUE)


pq <- ggplot(y,aes(x=dist/1e+3,y=qf)) +
  geom_pointrange(aes(ymin=lqf,ymax=uqf,colour="Prob")) +
  geom_pointrange(aes(x=(dist+10000)/1e+3,ymin=lqn*1e-3,ymax=uqn*1e-3,y=qn*1e-3,colour="Number")) +
  geom_path() + geom_hline(yintercept=1e-4,linetype="dashed",col="darkblue") +
  labs(x="Distance to TSS (kb)",y="Estimate of q") +
  geom_label_repel(aes(x=x,label=lab),
                   data=data.frame(x=1950,qf=1e-4,lab="q == 10^{-4}"),
                   hjust=1.2,vjust=-3,#fontface="bold",
                   col="darkblue",size=2.5,parse=TRUE) +
  facet_grid(. ~ what) +
  ## scale_y_log10() +
  scale_x_continuous(breaks=c(500,1000,1500,2000)) +
  scale_y_continuous(sec.axis = sec_axis(name = "Relative humidity [%]"))  +
  theme_pubr() +
  background_grid() +
  theme(axis.text=element_text(size=8))
pq

## annotate(geom="text",x=50,y=1e-4,label="q == 10^{-4}",vjust=-1,col="darkblue",parse=TRUE)

plot_grid(q,labels=c("a",""),ncol=1)
ggsave("~/fig-coloc-gtex.pdf",height=4,width=8/3)
plot_grid(p,labels=c("",""),ncol=1)
ggsave("~/fig-coloc-gtex-ncv.pdf",height=4,width=8/3)
