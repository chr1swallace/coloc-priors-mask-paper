#!/usr/bin/env Rscript

source("~/Projects/coloc-cond-mask/load-hapdata.R")

library(randomFunctions)
library(magrittr)
devtools::load_all("~/RP/coloc")
args <- getArgs(defaults=list(N=1000,NSIM=100,NCV=4,NSNP=250,SPECIAL=4,ld="highld",pop="eur"),
                numeric=c("N","NSIM","NCV"))

hdata <- hread(args$ld,args$pop)
LD <- hdata$LD
LD.alt <- hdata$LD.alt
dfsnps <- hdata$snps
h <- hdata$h
if(args$NSNP < nrow(dfsnps)) {
    dfsnps <- dfsnps[1:args$NSNP,]
    h <- h[,1:args$NSNP]
    LD <- LD[1:args$NSNP,1:args$NSNP]
}


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

R2 <- c(0.01,0.03,0.05)

addmethod <- function(L,meth) {
    dt <- sapply(L, is.data.table)  %>% which()
    if(!length(dt))
        return(NULL)
    dt <- rbindlist(L[dt])
    dt[,method:=meth]
    copy(dt)
}

plotter <- function(A,nCV) {
    p <- (A$beta/sqrt(A$varbeta))  %>% abs()  %>% pnorm(.,lower=FALSE)
    mycol <- ifelse(A$snp %in% nCV, "red","black")
    plot(seq_along(A$beta),-log10(p/2),col=mycol)
    ## points(nCV,-log10(p[nCV]/2),col="red")
}

    ## simulations options
## 5 = share 1, weakest effect for one, strongest for other
## 4 = share 2, opposite effects
## 3 = share 1, equal effect (weakest) + indep each
## 2 = share 2, equal effects
## 1 = share 1, equal effect (strongest), + indep each
## 0 = share 0
runone <- function(isim=0) {
    message("simulation ",isim)
    ## sample CVs. Ugly code, but works.
    if(args$NCV==4) { ## sample common CVs
        if(args$SPECIAL %in% c("2","4")) {
            CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),2)
            nCV1 <- nCV2 <- dfsnps$id[CV]
            beta1 <- beta2 <- sort(sample(c(1:9)/6,2))
            if(args$SPECIAL=="4")
                beta2 <- rev(beta2)
        } else if(args$SPECIAL=="0") {
            CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),4)
            nCV1 <- dfsnps$id[CV[1:2]]
            nCV2 <- dfsnps$id[CV[3:4]]
            beta1 <- sort(sample(c(1:9)/6,2))
            beta2 <- beta1 #sort(sample(c(1:6)/6,args$NCV))
        } else if(args$SPECIAL %in% c("1","3","5")) {
            CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),3)
            nCV1 <- dfsnps$id[CV[1:2]] # A B
            nCV2 <- dfsnps$id[CV[c(1,3)]] # A C
            b <- sort(sample(c(1:9)/6,2))
            if(args$SPECIAL=="1") {
                beta1 <- b[c(2,1)]
                beta2 <- b[c(2,1)]
            } else if(args$SPECIAL=="3") {
                beta1 <- b[c(1,2)]
                beta2 <- b[c(1,2)]
            } else if(args$SPECIAL=="5") {
                beta1 <- b[c(2,1)]
                beta2 <- b[c(1,2)]
            }
        } else {
            stop("special not programmed yet: ",args$SPECIAL)
        }
    } else if(args$NCV==3) {
        if(args$SPECIAL=="0") {
            CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),3)
            nCV1 <- dfsnps$id[CV[1:2]] # A, B
            nCV2 <- dfsnps$id[CV[3]] # C
            beta1 <- sort(sample(c(1:9)/6,2))
            beta2 <- beta1[2] #sort(sample(c(1:6)/6,1))
        } else {
            CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),2)
            nCV1 <- dfsnps$id[CV] # A, B
            nCV2 <- nCV1[1] # A
            b <- sort(sample(c(1:9)/6,2))
            if(args$SPECIAL=="1") {
                beta1 <- b[c(2,1)]
                beta2 <- b[c(2)] # share strong
            } else if(args$SPECIAL=="2") {
                beta1 <- b[c(1,2)]
                beta2 <- b[1] # share weak
            } else if(args$SPECIAL=="3") {
                beta1 <- b[c(2,1)]
                beta2 <- b[1] # share strong but with weak effect
            }
        }
    } else {
        stop("nCV not programmed yet: ",args$nCV)
    }


    ## CV=sample(which(dfsnps$AFR > 0.1 & dfsnps$AFR < 0.9),3)
    ##         nCV1 <- dfsnps$id[CV[1:2]] # A, B
    ##         nCV2 <- dfsnps$id[CV[3]] # C
    ##         beta1 <- sort(sample(c(1:9)/6,2))
    ##         beta2 <- beta1[2] #sort(sample(c(1:6)/6,1))
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
    usedata <- function(A) max(abs(A$beta)/sqrt(A$varbeta))  > 4.89
## stepwise regressions
    A1 <- makeD(y1,X1) # all signals, first dataset
    A2 <- makeD(y2,X2) # all signals, first dataset

    ## par(mfrow=c(2,1))
    ## plotter(A1,nCV1)
    ## plotter(A2,nCV2)
    
    
    ## indep analysis
    ## csignals1 <- finemap.indep.signals(A1,LD=LD,aligned=TRUE,zthr=4.89,maxhits=2)
    ## csignals2 <- finemap.indep.signals(A2,LD=LD,aligned=TRUE,zthr=4.89,maxhits=2)
    signals1 <- lapply(R2,function(r2) finemap.indep.signals(A1,LD=LD,zthr=4.89,r2thr=r2,maxhits=2))
    signals2 <- lapply(R2,function(r2) finemap.indep.signals(A2,LD=LD,zthr=4.89,r2thr=r2,maxhits=2))
    if(!length(signals1[[1]]) || !length(signals2[[1]]))
        return(NULL)
    
    dresults <- coloc.detail(A1,A2,p12=1e-6)  
      ## coloc.process(., hits1=signals1, LD=LD,r2thr=0.1^2)
      ## coloc.process(.,
      ##               hits1=names(signals1),
      ##               hits2=names(signals2),LD=LD)
    ## mask-mask with variable thresholds
    res <- lapply(R2, function(r2) { 
        coloc.signals(dataset1=c(A1,list(method="mask")),
                      dataset2=c(A2,list(method="mask")),LD=LD,
                      p12=1e-6,zthr=4.89,maxhits=2,r2thr=r2) %>%
          cbind(.,r2.process=r2,r2.finemap=r2)
    })  %>% addmethod("mask")
    ## res.alt <- lapply(R2, function(r2) { 
    ##         coloc.signals(c(A1,list(method="mask")),
    ##                       c(A2,list(method="mask")),LD=LD.alt,
    ##                       p12=1e-6,zthr=4.89,maxhits=2,r2thr=r2) %>%
    ##           cbind(.,r2.process=r2,r2.finemap=r2)
    ## })  %>% addmethod("mask.alt")
     ## result  <- lapply(R2, function(r2) {
    ##     lapply(seq_along(R2), function(i) {
    ##         if(all(signals1[[i]] <= 4.89) & all(signals2[[i]] <= 4.89))
    ##             next
    ##       coloc.process(dresults,
    ##                     hits1=names(signals1[[i]])[abs(signals1[[i]])>4.89],
    ##                     hits2=names(signals2[[i]])[abs(signals2[[i]])>4.89],
    ##                     LD=LD, r2thr=r2, p12=1e-6)  %>% cbind(.,r2.finemap=R2[[i]])
    ##     })  %>% rbindlist()  %>% cbind(.,r2.process=r2)
    ## })  %>% rbindlist()
    ## result[,method:="mask"]

    ## cond-mask, variable thresholds
    cres1 <- lapply(R2, function(r2) { 
            coloc.signals(c(A1,list(method="cond")),
                          c(A2,list(method="mask")),LD=LD,
                          p12=1e-6,zthr=4.89,maxhits=2,r2thr=r2) %>%
              cbind(.,r2.process=r2,r2.finemap=r2)
    })  %>% addmethod("condmask")
   ##  cres1.alt <- lapply(R2, function(r2) { 
   ##        coloc.signals(c(A1,list(method="cond")),
   ##                        c(A2,list(method="mask")),LD=LD.alt,
   ##                        p12=1e-6,zthr=4.89,maxhits=2,r2thr=r2) %>%
   ##            cbind(.,r2.process=r2,r2.finemap=r2)
   ##  }) %>% addmethod("condmask.alt")
     
   ## cond-mask, variable thresholds
    cres2 <- lapply(R2, function(r2) { 
            coloc.signals(c(A1,list(method="mask")),
                          c(A2,list(method="cond")),LD=LD,
                          p12=1e-6,zthr=4.89,maxhits=2,r2thr=r2) %>%
              cbind(.,r2.process=r2,r2.finemap=r2)
    })  %>% addmethod("maskcond")
    
   ##  cres2.alt <- lapply(R2, function(r2) { 
   ##          coloc.signals(c(A1,list(method="mask")),
   ##                        c(A2,list(method="cond")),LD=LD.alt,
   ##                        p12=1e-6,zthr=4.89,maxhits=2,r2thr=r2) %>%
   ##            cbind(.,r2.process=r2,r2.finemap=r2)
   ##  })  %>% addmethod("maskcond.alt")
    
     ## cond-cond, no thresholds
    ## cond.alt <- coloc.signals(A1,A2,LD=LD.alt,method="cond",p12=1e-6,zthr=4.89,maxhits=2)
    ## cond.alt$method <- "cond.alt"
    cond <- coloc.signals(A1,A2,LD=LD,method="cond",p12=1e-6,zthr=4.89,maxhits=2)
    cond$method <- "cond"
    ## cond1 <- coloc.signals(c(A1,list(method="cond")),A2,LD=LD,p12=1e-6,zthr=4.89,maxhits=2)
    ## cond1$method <- "cond1"
    ## cond2 <- coloc.signals(A1,c(A2,list(method="cond")),LD=LD,p12=1e-6,zthr=4.89,maxhits=2)
    ## cond2$method <- "cond2"
    ## cond1.alt <- coloc.signals(c(A1,list(method="cond")),A2,LD=LD.alt,p12=1e-6,zthr=4.89,maxhits=2)
    ## cond1.alt$method <- "cond1.alt"
    ## cond2.alt <- coloc.signals(A1,c(A2,list(method="cond")),LD=LD.alt,p12=1e-6,zthr=4.89,maxhits=2)
    ## cond2.alt$method <- "cond2.alt"

    ## single
    single <- coloc.process(dresults)  %>% as.data.table()
    single[,hit1:=names(signals1[[1]])[[1]]]
    single[,hit2:=names(signals2[[1]])[[1]]]
    single[,method:="single"]
    z1 <- structure(A1$beta/sqrt(A1$varbeta),names=A1$snp)
    z2 <- structure(A2$beta/sqrt(A2$varbeta),names=A2$snp)
    single[,hit1.margz:=z1[hit1]]
    single[,hit2.margz:=z2[hit2]]

    result <- rbind(single,
                    res,cres1,cres2,
                    cond,#cond1,cond2,
                    ## res.alt,cres1.alt,cres2.alt,
                    ## cond.alt,cond1.alt,cond2.alt,
                    fill=TRUE)
    result[,r2.cvA1.hits:=LD[ hit1,nCV1[1] ]^2 ]
    result[,r2.cvB1.hits:=LD[ hit1,nCV1[2] ]^2 ]
    result[,r2.cvA2.hits:=LD[ hit2,nCV2[1] ]^2 ]
    if(args$NCV==4)
        result[,r2.cvB2.hits:=LD[ hit2,nCV2[2] ]^2 ]
    
    ## conditional analysis
    result$special <- args$SPECIAL
    result$r2.tr1 <- if(length(nCV1)==1) { NA } else { LD[nCV1,nCV1][1,2]^2 }
    result$r2.tr2 <- if(length(nCV2)==1) { NA } else { LD[nCV2,nCV2][1,2]^2 }
    result$r2.tr12 <- max(LD[nCV1,nCV2]^2)
    result[,ntested:=.N,by=c("method","r2.finemap","r2.process")]
    
    ## result$details <- paste(c(nCV1,beta1,nCV2,beta2#,ld[1,2],ld[3,4],max(ld[1:2,3:4])
    ##                           ),collapse="/")
    result$isim <- isim
    result
}

library(parallel)
options(mc.cores=5)
results <- mclapply(1:args$NSIM,runone)
nulls <- sapply(results,is.null)
results  <- rbindlist(results[!nulls])
results$NCV <- if(args$NCV==3) { "onetwo" } else { "twotwo" }
results$ld <- args$ld
results[,N:=args$N]
results$NSNP <- args$NSNP


patt <- "cvaryr2-v13"
tmp <- tempfile(tmpdir=d,pattern=patt,fileext=".RData")
save(results,file=tmp)

message("COMPLETE")
