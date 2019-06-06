#!/usr/bin/env Rscript
library(randomFunctions)
library(snpStats)
library(magrittr)
library(data.table)
source("~/DIRS.txt")

## input data
dref <- file.path(REFDATA,"1000GP_Phase3")
## store output
d <- COLOCINDEP
dhg <- file.path(d,"input")



if(FALSE) { # generate files
    
samples <- fread(file.path(dref,"1000GP_Phase3.sample"))
leg <- fread(cmd=paste("zcat",
                   file.path(dref,"1000GP_Phase3_chr10.legend.gz")))
## https://stackoverflow.com/questions/16911773/collapse-runs-of-consecutive-numbers-to-ranges
findIntRuns <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x)==1) as.character(x) else paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

## what we need to specify
## high ld - il2ra region
## 10p-6030000-6220000
    whs <- which(samples$GROUP=="AFR")
    whs <- c(whs*2,whs*2-1)  %>% sort()
cuts.afr <- findIntRuns(whs)  %>% paste(.,collapse=",")
    whs <- which(samples$GROUP=="EUR")
    whs <- c(whs*2,whs*2-1)  %>% sort()
cuts.eur <- findIntRuns(whs)  %>% paste(.,collapse=",")
mylast <- function(x) { x[length(x)]  }
minl <- which(leg$position > 6030000)[1]
maxl <- which(leg$position < 6220000)  %>% mylast()
    whl <- which(leg$position > 6030000 &
                  leg$position < 6220000 &
                  pmin(leg$EUR,1-leg$EUR) > 0.01) # 1% MAF
    ## cutl <- findIntRuns(whl)  %>% paste0(" -e ",.,"p")  %>% paste(.,collapse="")
    
## rows = variants, no header
## cols = samples, first col is snp label
system(paste0("zcat ",
            file.path(dref,"chr10.hap.gz"),
            " | sed -n ",minl,",",maxl,"p",
            " | cut -d' ' -f ",cuts.eur,
            "> ", file.path(dhg,"highld-eur.hap")))
system(paste0("zcat ",
            file.path(dref,"chr10.hap.gz"),
            " | sed -n ",minl,",",maxl,"p",
            " | cut -d' ' -f ",cuts.afr,
            "> ", file.path(dhg,"highld-afr.hap")))
## rows = variants, 1st row is header
com <- paste0("zcat ",
             file.path(dref,"1000GP_Phase3_chr10.legend.gz"),
             "| sed -n -e 1p -e ",minl+1,",",maxl+1,"p ",
             "> ", file.path(dhg,"highld.leg"))
com
system(com)

## what we need to specify
## low ld - 1st 5000 snps
whs <- which(samples$GROUP=="AFR")
    whs <- c(whs*2,whs*2-1)  %>% sort()
cuts.afr <- findIntRuns(whs)  %>% paste(.,collapse=",")
    whs <- which(samples$GROUP=="EUR")
    whs <- c(whs*2,whs*2-1)  %>% sort()
cuts.eur <- findIntRuns(whs)  %>% paste(.,collapse=",")
minl <- 1 #which(leg$position > 6030000)[1]
maxl <- 12000 #which(leg$position < 6220000)  %>% mylast()

## rows = variants, no header
## cols = samples, first col is snp label
system(paste0("zcat ",
            file.path(dref,"chr10.hap.gz"),
            " | sed -n ",minl,",",maxl,"p",
            " | cut -d' ' -f ",cuts.eur,
            "> ", file.path(dhg,"lowld-eur.hap")))
system(paste0("zcat ",
            file.path(dref,"chr10.hap.gz"),
            " | sed -n ",minl,",",maxl,"p",
            " | cut -d' ' -f ",cuts.afr,
            "> ", file.path(dhg,"lowld-afr.hap")))
## rows = variants, 1st row is header
com <- paste0("zcat ",
             file.path(dref,"1000GP_Phase3_chr10.legend.gz"),
             "| sed -n -e 1p -e ",minl+1,",",maxl+1,"p ",
             "> ", file.path(dhg,"lowld.leg"))
system(com)

}

################################################################################

## functions to make h
hread <- function(ld=c("lowld","highld"),pop=c("eur","afr")) {
    ld <- match.arg(ld)
    snps <- fread(paste0(dhg,"/",ld,".leg"))
    ## map <- fread(file.path(d,"example/ex.map"))
### h <- fread(file.path(d,"example/ex.haps"))[1:1000,]
    ## snps <- fread(file.path(d,"example/ex.leg"))
## map <- fread(file.path(d,"example/ex.map"))
    ## wh <- which(c(0,diff(map[["Genetic_Map(cM)"]]))>0.1)
    f <- function(h) {
        h <- as.matrix(h)
        h1 <- h[, seq(1,ncol(h)-1,by=2)] #[,samples$GROUP=="AFR"]
        h2 <- h[, seq(2,ncol(h),by=2)] #[,samples$GROUP=="AFR"]
        t(cbind(h1,h2))
    }
    h.eur <- fread(paste0(dhg,"/",ld,"-eur.hap"))  %>% f()
    h.afr <- fread(paste0(dhg,"/",ld,"-afr.hap"))  %>% f()
    h <- if(pop=="eur") { h.eur } else { h.afr }
    h.alt <-  if(pop=="eur") { h.afr } else { h.eur }
    maf  <-  colMeans(h)
    maf.alt  <-  colMeans(h.alt)
    use <- maf > 0.01 & maf < 0.99 & maf.alt > 0.01 & maf.alt < 0.99 & apply(h,2,var)>0
    h <- h[,use,drop=FALSE] 
    h.alt <- h.alt[,use,drop=FALSE] 
    dfsnps <- snps[use,,drop=FALSE]
    dfsnps$id <- make.names(dfsnps$id)
    dfsnps$maf <- colMeans(h)
    dfsnps$maf.alt <- colMeans(h.alt)
    LD <- cor(h)
    LD.alt <- cor(h.alt)
    ## LD <- as.matrix(make.positive.definite(LD))
    dimnames(LD) <- dimnames(LD.alt)  <- list(dfsnps$id,dfsnps$id)
    return(list(h=h,snps=dfsnps,LD=LD,LD.alt=LD.alt))
}
