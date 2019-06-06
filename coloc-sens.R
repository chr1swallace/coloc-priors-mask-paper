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

addmethod <- function(L,meth) {
    dt <- sapply(L, is.data.table)  %>% which()
    if(!length(dt))
        return(NULL)
    dt <- rbindlist(L[dt])
    dt[,method:=meth]
    copy(dt)
}

plotter <- function(obj,wh,position=1:nrow(obj$df)) {
    znm <- if(wh==1) { "z.df1" } else {"z.df2" }
    p <- pnorm(abs(obj$df[[znm]]),lower=FALSE)*2
    ## mycol <- ifelse(A$snp %in% nCV, "red","black")
    Pal <- colorRampPalette(c('white','blue'))

#This adds a column of color values
# based on the y values
    Col <- Pal(100)[ceiling(100*obj$df$H4)]
    plot(position,-log10(p),col="gray20",
           bg = Col, # Fill colour
                pch = 21, # Shape: circles that can filed
                frame.plot = FALSE, # Remove the frame 
               xlab="Chromosome position",ylab="-log10(p)",xaxt='n')
   axis(side=1,labels=FALSE) 
  ## points(nCV,-log10(p[nCV]/2),col="red")
}

    ## simulations options
## 5 = share 1, weakest effect for one, strongest for other
## 4 = share 2, opposite effects
## 3 = share 1, equal effect (weakest) + indep each
## 2 = share 2, equal effects
## 1 = share 1, equal effect (strongest), + indep each
## 0 = share 0
simone <- function() {
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
    list(A1,A2)
}



set.seed(43)
A <- simone()
dresults <- coloc.detail(A[[1]],A[[2]],p12=1e-6)  
m <- match(dresults$df$snp, dfsnps$id)
pdf("~/fig-coloc-sens-1.pdf",height=6,width=8)
layout(mat = matrix(1:4,2,2),
       heights = c(1, 1), # Heights of the two rows
       widths = c(2, 3)) # Widths of the two columns
par(mar = c(3, 3, 2, 1), # Dist' from plot to side of page
    mgp = c(2, 0.4, 0), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.01) # Reduce tick length
plotter(dresults,1,dfsnps$position[m]); title(main="Manhattan plot, trait 1",adj=0)
plotter(dresults,2,dfsnps$position[m]); title(main="Manhattan plot, trait 2",adj=0)
sensitivity.coloc(dresults,preserve.par=TRUE)
dev.off()

## only top end of p12 passes
set.seed(403)
A <- simone()
dresults <- coloc.detail(A[[1]],A[[2]],p12=1e-6)  
pdf("~/fig-coloc-sens-2.pdf",height=6,width=8)
layout(mat = matrix(1:4,2,2),
       heights = c(1, 1), # Heights of the two rows
       widths = c(2, 3)) # Widths of the two columns
    par(mar = c(3, 3, 2, 1), # Dist' from plot to side of page
    mgp = c(2, 0.4, 0), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.01 # Reduce tick length
    )
plotter(dresults,1,dfsnps$position[m]); title(main="Manhattan plot, trait 1",adj=0)
plotter(dresults,2,dfsnps$position[m]); title(main="Manhattan plot, trait 2",adj=0)
sensitivity.coloc(dresults,preserve.par=TRUE)
dev.off()
