#!/usr/bin/env Rscript
library(randomFunctions)
library(snpStats)
library(magrittr)
devtools::load_all("~/RP/coloc")
source("~/DIRS.txt")
args <- getArgs(defaults=list(N=1000,NSIM=100,NCV=2,SPECIAL=2),numeric=c("N","NSIM","NCV"))
library(data.table)

d <- COLOCINDEP

## FIRST 5000 snps in Africans only
dref <- file.path(REFDATA,"1000GP_Phase3")
fhg <- file.path(d,"input")#tempfile(tmpdir=dtmp)
## h <- fread(paste("zcat",
##             file.path(dref,"chr21_afr.hap.gz"),
##             "|head -n 100"))
## crange <- paste0("7-",ncol(h)+6)
## system(paste("zcat",
##             file.path(dref,"chr21_afr.hap.gz"),
##             "|head -n 5000 | cut -d' ' -f",crange,
##             "> ", paste0(fhg,"/hap")))
## system(paste("zcat",
##              file.path(dref,"1000GP_Phase3_chr21.legend.gz"),
##              "| head -n 5001 > ", paste0(fhg,"/leg")))

library(corpcor)
## read in haplotypes
h <- fread(paste0(fhg,"/hap"))
h <- as.matrix(h)
snps <- fread(paste0(fhg,"/leg"))
## samples <- fread(file.path(dref,"1000GP_Phase3.sample"))
## map <- fread(file.path(d,"example/ex.map"))
## wh <- which(c(0,diff(map[["Genetic_Map(cM)"]]))>0.1)

h1 <- h[, seq(1,ncol(h)-1,by=2)] #[,samples$GROUP=="AFR"]
h2 <- h[, seq(2,ncol(h),by=2)] #[,samples$GROUP=="AFR"]
h <- t(cbind(h1,h2))
use <- snps$AFR>0.01 & snps$AFR < 0.99 & apply(h,2,var)>0
h <- h[,use,drop=FALSE]
dfsnps <- snps[use,,drop=FALSE]
dfsnps$maf <- colMeans(h)
dfsnps$id <- make.names(dfsnps$id)
LD <- cor(h)
## LD <- as.matrix(make.positive.definite(LD))
dimnames(LD) <- list(dfsnps$id,dfsnps$id)

## XX <- new("SnpMatrix", as.matrix(freq))S""

makeD <- function(y,X,best=NULL) {
    if(is.null(best)) {
        m <- snp.rhs.estimates(y~1,snp.data=X,family="Gaussian")
    } else {
        df <- data.frame(y=y,best=as(X[,best],"numeric"),row.names=rownames(X))
        m <- snp.rhs.estimates(y ~ .,data=df,snp.data=X,family="Gaussian")
    }
    nulls <- sapply(m,is.null)
    if(any(nulls)) 
        m <- m[!nulls]
    b <- sapply(m,"[[","beta")
    v <- sapply(m,"[[","Var.beta")
    z <- b/sqrt(v)
    list(N=args$N,
         MAF=dfsnps$maf[!nulls],
         beta=b,
         varbeta=v,
         type="quant",
         sdY=sd(y),
         snp=dfsnps$id[!nulls])
}
getmaxz <- function(D) {
    wh <- which.max(abs(D$beta)/sqrt(D$varbeta))
    D$snp[ wh ]
}
valmaxz <- function(D) {
    z <- abs(D$beta)/sqrt(D$varbeta)
    max(z)
}

## simulations options
## 3 = share 1, stronger effect
## 2 = share 1, weaker effect
## 1 = share 1, equal effect (strongest), + indep each
## 0 = share 0
runone <- function(isim=0) {
    message("simulation ",isim)
    ## sample common CVs
    if(args$SPECIAL=="0") {
        CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),3)
        nCV1 <- dfsnps$id[CV[1:2]]
        nCV2 <- dfsnps$id[CV[3]]
        beta1 <- sort(sample(c(1:6)/6,2))
        beta2 <- sort(sample(c(1:6)/6,1))
    } else {
        CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),2)
        nCV1 <- dfsnps$id[CV]
        nCV2 <- nCV1[1]
        b <- sort(sample(c(1:6)/6,3))
        if(args$SPECIAL=="1") {
            beta1 <- b[c(2,3)]
            beta2 <- b[c(3)]
        } else if(args$SPECIAL=="2") {
            beta1 <- b[c(1,2)]
            beta2 <- b[3]
        } else if(args$SPECIAL=="3") {
            beta1 <- b[c(2,3)]
            beta2 <- b[1]
        }
    }    

    ## make genotypes
    nr <- nrow(h)
    G1 <- h[sample(1:nr,args$N,replace=TRUE),] + h[sample(1:nr,args$N,replace=TRUE),]
    G2 <- h[sample(1:nr,args$N,replace=TRUE),] + h[sample(1:nr,args$N,replace=TRUE),]
    colnames(G1) <- colnames(G2) <- dfsnps$id
    rownames(G1) <- rownames(G2) <- paste0("I",1:args$N)
    X1 <- new("SnpMatrix",G1+1)
    X2 <- new("SnpMatrix",G2+1)

    ## outcomes
    y1  <-  rnorm(args$N) + tcrossprod(G1[,nCV1],t(beta1))
    y2  <-  rnorm(args$N) + tcrossprod(G2[,nCV2],t(beta2))
#hist(y,breaks=100)

## stepwise regressions
    A1 <- makeD(y1,X1) # all signals, first dataset
    bestA1 <- getmaxz(A1)
    B1 <- makeD(y1,X1,bestA1) # second signal, first dataset
    bestB1 <- getmaxz(B1)
    C1 <- makeD(y1,X1,bestB1) # first signal, first dataset
    A2 <- makeD(y2,X2) # all signals, second dataset
                 
    ## indep analysis
    signals1 <- finemap.indep.signals(A1,LD=LD,zthr=qnorm(5e-2/2,lower.tail=FALSE))[1:length(beta1)] 
    signals2 <- finemap.indep.signals(A2,LD=LD,zthr=qnorm(5e-2/2,lower.tail=FALSE))[1:length(beta2)]
    print(signals1)
    print(signals2)
    mask <- coloc.detail(A1,A2)  %>%
      coloc.process(., hits1=signals1, LD=LD)
      ## coloc.process(., hits1=signals1, hits2=signals2, LD=LD)

    ## stepwise analysis - conditioning each time
    cond1 <- coloc.detail(C1,A2)  %>% coloc.process()
    cond2 <- coloc.detail(B1,A2)  %>%  coloc.process()
    ## stepwise analysis - conditioning second time
    step1 <- coloc.detail(A1,A2)  %>% coloc.process()
    step2 <- cond2 #coloc.detail(B1,A2)  %>%  coloc.process()
    ## vmz <- c(valmaxz(B1))
    step <- rbind(step1,step2)
    cond <- rbind(cond1,cond2)
                  ## if(vmz[1] >4) { stepB1 } else { NULL })
    
    mask$method <- "mask"
    step$method <- "step"
    cond$method <- "cond"
    step$mask <- c("",bestA1)
    cond$mask <- c(bestB1,bestA1)
    mask$mask <- rev(signals1)
    cond$r2.cvA.hits <- step$r2.cvA.hits <- LD[nCV1[1],c(bestA1,bestB1)]
    cond$r2.cvB.hits <- step$r2.cvB.hits <- LD[nCV1[2],c(bestA1,bestB1)]
    mask$r2.cvA.hits=LD[nCV1[1],signals1]
    mask$r2.cvB.hits=LD[nCV1[2],signals1]
    
    result <- as.data.table(rbind(mask,step,cond))
    ## result[,best1:=as.character(best1)]
    ## result[,best2:=as.character(best2)]
    ## result[,r2.tr1.cv1:=LD[nCV1,nCV1]^2[1,2]]
    ## result[,r2.tr1.cv2:=LD[best1,nCV1[2]]^2]
    ## result[,r2.tr2.cv:=LD[best1,nCV2]^2]

    result$r2.tr1 <- LD[nCV1,nCV1][1,2]^2
    result$special <- args$SPECIAL
    result$details <- paste(c(nCV1,beta1,nCV2,beta2),collapse="/")
    result$isim <- isim
    result
}

results <- lapply(1:args$NSIM,runone)  %>% rbindlist()
## results <- lapply(1:2,runone)  %>% rbindlist()
    
tmp <- tempfile(tmpdir=d,pattern="c12sim",fileext=".RData")
save(results,file=tmp)
