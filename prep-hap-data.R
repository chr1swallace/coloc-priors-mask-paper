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

if(FALSE) {
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
  maxl <- which(leg$position > 6030000)[2000]  %>% mylast()
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
​
​
​
######################
​
## 2. Convert these to lists of h, snps and LD
​
hread <- function(ld=c("low","medium","medium2","high")) {
  if(ld == "low"){
    ld.int = "highld"
    pop = "afr"
  } else {
    if(ld == "medium"){
      ld.int = "highld"
      pop = "eur"
    } else {
      if(ld == "medium2"){
        ld.int = "lowld"
        pop = "afr"
      } else {
        ld.int = "lowld"
        pop = "eur"        
      }
    }
  }
  snps <- fread(paste0(dhg,"/",ld.int,".leg"))
  ## map <- fread(file.path(d,"example/ex.map"))
  ### h <- fread(file.path(d,"example/ex.haps"))[1:1000,]
  ## snps <- fread(file.path(d,"example/ex.leg"))
  ## map <- fread(file.path(d,"example/ex.map"))
  ## wh <- which(c(0,diff(map[["Genetic_Map(cM)"]]))>0.1)
  f <- function(h) {
    h <- as.matrix(h)
    h1 <- h[, seq(1,ncol(h)-1,by=2)] #[,samples$GROUP=="AFR"] ODD
    h2 <- h[, seq(2,ncol(h),by=2)] #[,samples$GROUP=="AFR"] EVEN
    t(cbind(h1,h2))
  }
  h.eur <- fread(paste0(dhg,"/",ld.int,"-eur.hap"))  %>% f()
  h.afr <- fread(paste0(dhg,"/",ld.int,"-afr.hap"))  %>% f()
  h <- if(pop=="eur") { h.eur } else { h.afr }
  maf  <-  colMeans(h)
  use <- maf > 0.01 & maf < 0.99 & apply(h,2,var)>0
  h <- h[,use,drop=FALSE] 
  dfsnps <- snps[use,,drop=FALSE]
  dfsnps$id <- make.names(dfsnps$id)
  dfsnps$maf <- colMeans(h)
  LD <- cor(h)
  ## LD <- as.matrix(make.positive.definite(LD))
  dimnames(LD) <- list(dfsnps$id,dfsnps$id)
  return(list(h=h,snps=dfsnps,LD=LD))
}
​
tmp <- hread(ld ="medium")
​
library(corrplot)
corrplot(tmp$LD^2, "color", "upper", tl.pos = "n")
​
######################
​
## 3. Make RDS files of h for code
​
low <- hread(ld = "low")
medium <- hread(ld = "medium")
medium2 <- hread(ld = "medium2")
high <- hread(ld = "high")

par(mfrow=c(2,2))
corrplot(low$LD^2, "color","upper", tl.pos = "n", cl.pos="n")
corrplot(medium$LD^2, "color","upper", tl.pos = "n", cl.pos="n")
corrplot(medium2$LD^2, "color","upper", tl.pos = "n", cl.pos="n")
corrplot(high$LD^2, "color","upper", tl.pos = "n", cl.pos="n")

summary(as.vector(low$LD^2))
summary(as.vector(medium$LD^2))
summary(as.vector(medium2$LD^2))
summary(as.vector(high$LD^2))
​
saveRDS(medium2, "~/scratch/lowld.RDS")
saveRDS(high, "~/scratch/highld.RDS")
