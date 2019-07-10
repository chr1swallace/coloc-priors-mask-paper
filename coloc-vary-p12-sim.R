#!/usr/bin/env Rscript

source("~/Projects/coloc-cond-mask/load-hapdata.R")

library(randomFunctions)
library(magrittr)
devtools::load_all("~/RP/coloc")
args <- getArgs(defaults=list(N=1000,NSNP=500,NSIM=100,NCV=3,SPECIAL=2,ld="highld",pop="eur"),
                numeric=c("N","NSIM","NCV","NSNP"))

hdata <- hread(args$ld,args$pop)
LD <- hdata$LD
dfsnps <- hdata$snps
h <- hdata$h
if(args$NSNP < nrow(dfsnps)) {
    dfsnps <- dfsnps[1:args$NSNP,]
    h <- h[,1:args$NSNP]
    LD <- LD[1:args$NSNP,1:args$NSNP]
}


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
    z <- abs(D$beta)/sqrt(D$varbeta)
    wh <- which.max(z)
    structure(z[wh],names=D$snp[ wh ])
}
valmaxz <- function(D) {
    z <- abs(D$beta)/sqrt(D$varbeta)
    max(z)
}

## simulations options
# args$SPECIAL=H
runone <- function(isim=0) {
    message("simulation ",isim)
    ## sample CVs. Ugly code, but works.
    pr.var1 <- rnorm(1,0,0.15)  %>% abs()
    pr.var2 <- rnorm(1,0,0.15)  %>% abs()
        
   if(args$SPECIAL %in% c("0","1","2","3")) { # share 0
       CV=sample(which(dfsnps$maf > 0.1 & dfsnps$maf < 0.9),2)
       nCV1 <- dfsnps$id[CV[1]]
       nCV2 <- dfsnps$id[CV[2]]
       beta1 <- beta2 <- sample(c(1:9)/6,1)
   } else if(args$SPECIAL=="4") {
       CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$maf < 0.9),1)
       nCV1 <- nCV2  <- dfsnps$id[CV[1]]
       beta1 <- beta2 <- sample(c(1:9)/6,1)
   } else {
       stop("special not programmed yet: ",args$SPECIAL)
   }
    if(args$SPECIAL %in% c("1","0")) {
        beta2 <- 0
        pr.var2 <- 0
    }
    if(args$SPECIAL %in% c("2","0")) {
        beta1 <- 0
        pr.var1 <- 0
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
    if(pr.var1>0) {
        b1 <- tcrossprod(G1[,nCV1],t(beta1))
        y1 <- sqrt(pr.var1) * b1/sd(b1) + rnorm(args$N,sd=sqrt(1-pr.var1))
    } else {
        y1 <- rnorm(args$N,sd=1)
    }
    
    if(pr.var2>0) {
        b2 <- tcrossprod(G2[,nCV2],t(beta2))
        y2 <- sqrt(pr.var2) * b2/sd(b2) + rnorm(args$N,sd=sqrt(1-pr.var2))
    } else {
        y2 <- rnorm(args$N,sd=1)
    }
    
#hist(y,breaks=100)

    usedata <- function(A) max(abs(A$beta)/sqrt(A$varbeta))  > 4.89
## stepwise regressions
    ## print(summary(y1))
    A1 <- makeD(y1,X1) # all signals, first dataset
    bestA1 <- getmaxz(A1)
                 
    ## print(summary(y2))
    A2 <- makeD(y2,X2) # all signals, first dataset
    bestA2 <- getmaxz(A2)
                 
    ## stepwise analysis - conditioning second time
    cd <- coloc.detail(A1,A2,p12=1e-4)
    ## prior.sens(cd)

    ## coloc.detail(A1,A2,p12=1e-4)$summary
    ## cd$summary
    ## prior.sens(cd, "H4>H3")

    pp <- 1e-6
    result <- #lapply(c(1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8), function(pp) {
        cbind(data.table(p12=pp),coloc.process(cd,p12=pp,p1=1e-4-pp,p2=1e-4-pp))
    ## })  %>% do.call("rbind",.)
    result[,hit1:=NULL][,hit2:=NULL]
    result[,best1:=unlist(best1)]
    result[,best2:=unlist(best2)]
    ## result$hit1 <- names(bestA1)
    ## result$hit2 <- names(bestA2)
    result$z1 <- bestA1
    result$z2 <- bestA2
                  ## if(vmz[1] >4) { stepB1 } else { NULL })
    ## cond$r2.cvA1.hits <- step$r2.cvA1.hits <- LD[nCV1[1],c(bestA1,bestB1)]
    ## cond$r2.cvB1.hits <- step$r2.cvB1.hits <- LD[nCV1[2],c(bestA1,bestB1)]
    ## cond$r2.cvA2.hits <- step$r2.cvA2.hits <- LD[nCV2[1],c(bestA2,bestA2,bestB2,bestB2)]
    ## cond$r2.cvB2.hits <- step$r2.cvB2.hits <- LD[nCV2[2],c(bestA2,bestA2,bestB2,bestB2)]
    ## mask$r2.cvA1.hits=LD[nCV1[1],signals1]
    ## mask$r2.cvB1.hits=LD[nCV1[2],signals1]
    ## mask$r2.cvA2.hits=LD[nCV2[1],rep(signals2,each=2)]
    ## mask$r2.cvB2.hits=LD[nCV2[2],rep(signals2,each=2)]

    result[,r2.cvA1.hits:=LD[ best1,nCV1[1] ]^2 ]
    result[,r2.cvA2.hits:=LD[ best2,nCV2[1] ]^2 ]
    
 ## stepwise analysis
    result$r2.tr12 <- max(LD[nCV1,nCV2]^2)
    
    ## result$details <- paste(c(nCV1,beta1,nCV2,beta2#,ld[1,2],ld[3,4],max(ld[1:2,3:4])
    ##                           ),collapse="/")
    result$isim <- isim
    result
}


results <- lapply(1:args$NSIM,runone)  %>% rbindlist()
results[,N:=args$N]
results$special <- args$SPECIAL
results$ld <- args$ld
results$pop <- args$pop
results$NSNP <- args$NSNP

patt <- "cvaryp12-v4"
tmp <- tempfile(tmpdir=d,pattern=patt,fileext=".RData")
save(results,file=tmp)
