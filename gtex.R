#!/usr/bin/env Rscript

d <- "/home/cew54/share/Data/gtex"
library(magrittr)
library(randomFunctions)
args <- getArgs(defaults=list(chr=14),
                numeric="taskid")
## of <- file.path("~/scratch/gtex-gfm",paste0(paste(args$chr,args$taskid,sep="-"),".rds"))
## if(file.exists(of))
##     stop("output file already exists: ",of)

SIZE=100
message("\n! ---------------------------------------")
message("! ",date())
message("! gtex-eqtl.R running with args:")
print(args)
message("! ---------------------------------------\n")


PLINK="/home/cew54/localc/bin/plink" # plink binary
BCFTOOLS="/home/cew54/localc/bin/bcftools" # bcftools binary

d <- "/home/cew54/share/Data/gtex/eur-Whole_Blood-001"

## load gt data
library(annotSnpStats)
(load(file.path(d, paste0("chr",args$chr,".RData"))))
dim(X)
head(X@samples)
head(X@snps)

## covariates
library(data.table)
covar <- fread(file.path(d, "../GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"))
m <- as.matrix(covar[,-1])
rownames(m) <- covar$ID
m <- t(m)
m <- m[rownames(X),]
m <- as.data.frame(m)

## expression
library(data.table)
expr <- fread(paste0("zcat ",file.path(d, "../GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz")))
## expr <- expr[,c("gene_id",rownames(X))]
setnames(expr,"#chr","chr")
expr <- expr[chr==paste0("chr",args$chr),]
dim(expr)
E <- as.matrix(expr[,-c(1:4)])
rownames(E) <- expr$gene_id
E <- t(E)[rownames(X),]
expr <- expr[,1:4]
## limit snps

cs0 <- col.summary(X)
wh <- which(is.na(cs0[,"z.HWE"]) |
            cs0[,"MAF"]<0.01 |
            cs0[,"Call.rate"]<0.99 |
            cs0[,"Certain.calls"]<0.75 | 
            abs(cs0[,"z.HWE"])>4)
if(length(wh)) {
    message("Dropping ",length(wh)," SNPs with |z.HWE|>5, MAF < 0.001 or call rate <0.99")
    X <- X[,-wh]
}
dim(X)

library(GUESSFM)
W <- 2e+6
TAG.R2=0.95
NSWEEP=500000
NSAVE=5000
NCHAINS=5
NEXP=3
COMFILE <- paste0("/home/cew54/scratch/gtex-gfm/chr",args$chr,"/runme.sh")
#if(file.exists(COMFILE))
#    unlink(COMFILE)

message("Expecting jobs: ",ncol(E))
files=list.files( paste0("/home/cew54/scratch/gtex-gfm/chr",args$chr))
nfiles=length(files)
if(file.exists(COMFILE)) {
comms=scan(COMFILE,what="",sep="\n")
ncom=length(comms)
nfiles=nfiles - 1
} else {
ncom=0
}
message("Job dirs found: ",nfiles)
message("Job commands found: ",ncom)

if(ncom==ncol(E) & !interactive())
    q("no")
    
for(j in 1:ncol(E)) {
    cat(j,"\r")
    outd <- paste0("/home/cew54/scratch/gtex-gfm/chr",args$chr,"/",j)
    if(!file.exists(outd))
        dir.create(outd,recursive = TRUE)
    
    f.data <- file.path(outd,"data.RData")
    f.par <- file.path(outd, "par.xml")
    if(file.exists(f.par) ) {
#        message("output file already exists, skipping: ",f.par)
	next()
    }
    use <- which( (X@snps$position > expr$start[j] - W) & (X@snps$position < expr$start[j] + W) )
    sX <- sm(X)[,use]
    tags <- tag(sX, method="single", tag.threshold=TAG.R2)
    message("after tagging at ",TAG.R2,", matrix now has ",length(unique(tags(tags)))," SNPs.")
    ## save(tags, file=f.tags)
    DATA <- sX[, unique(tags(tags))]
    ## save data
    y <- E[,j]
    snps <- X@snps[colnames(sX),]
    thisgene <- expr[j,]
    save(thisgene,snps,file=f.data)
    message("Samples: ", length(y))
    message("SNPs: ",ncol(DATA))
    coms <- run.bvs(X=DATA,Y=y, covars=m,
                    gdir=outd,nsweep=NSWEEP, family="gaussian",tag.r2=NA,
                    nsave=NSAVE,nchains=NCHAINS,nexp=NEXP,run=FALSE)
    cat(coms,file=COMFILE,sep="\n",append=TRUE)
}

