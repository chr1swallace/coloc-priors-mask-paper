
library(data.table)

setwd("~/scratch/gtex")
## system("wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz")
## system("tar zxvf GTEx_Analysis_v7_eQTL.tar.gz")

## signif signals
x=fread(cmd="zcat GTEx_Analysis_v7_eQTL/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz")
x[,c("chr","pos.37","a1","a2","b37"):=tstrsplit(variant_id,"_")]
x[,pos.37:=as.numeric(pos.37)]

## ## liftover
## library(rtracklayer)
## library(GenomicRanges)
## library(magrittr)
## path = file.path("~/share/Data/reference/hg19ToHg38.over.chain")
## system(paste0("gunzip ",path,".gz"))
## ch = import.chain(path)
## ch
## system(paste0("gzip ",path))
## tmp <- GRanges(seqnames=Rle(paste0("chr",x$chr)),
##                ranges=IRanges(start=x$pos.37,width=1))
## tmp38 = (liftOver(tmp, ch))
## class(tmp38)

## ln <- lengths(tmp38)
## table(ln)
## tmp38 <- tmp38[which(ln==1)]  %>% unlist()
## x <- x[-which(ln!=1),]
## x[,pos.38 := start(tmp38)]

## x[,id38:=paste(chr,pos.38,a1,a2,sep=":")]
## cat(x$id38,file="gtex-snps.txt",sep="\n")

## LD - not this - look at simulations
library(snpStats)
library(magrittr)
d <- "/home/cew54/share/Data/reference/1000GP_Phase3"

## EUR samples only
## https://stackoverflow.com/questions/16911773/collapse-runs-of-consecutive-numbers-to-ranges
source("~/DIRS.txt")
dref <- file.path(REFDATA,"1000GP_Phase3")
samples <- fread(file.path(dref,"1000GP_Phase3.sample"))
findIntRuns <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x)==1) as.character(x) else paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}
  whs <- which(samples$GROUP=="EUR")
    whs <- c(whs*2,whs*2-1)  %>% sort()
cuts.eur <- findIntRuns(whs+1)  %>% paste(.,collapse=",") # +1 because we will add a line number

devtools::load_all("~/RP/coloc")
windows <- seq(1e+5,1e+6,by=1e+5)

## chr 22
## RES <- lapply(1:22, function(ichr) {
for(ichr in 1:22) {

    outfile <- paste0("map3.",ichr)
    if(file.exists(outfile))
        next
    
    message("chromosome ",ichr)
    xx <- x[chr==ichr,]
    xx[,id:=paste(chr,pos.37,a1,a2,sep=":")]
    legfile <- file.path(dref,paste0("1000GP_Phase3_chr",ichr,".legend.gz"))
    hapfile <- file.path(dref,paste0("1000GP_Phase3_chr",ichr,".hap.gz"))

    message("reading legend ",basename(legfile))
    leg <- fread(cmd=paste("zcat",legfile))
    leg[,id37:=paste(ichr,position,a0,a1,sep=":")]
    ## table(xx$pos.37 %in% leg$position)
    table(xx$id %in% leg$id37) # bad - need liftover
    m <- match(xx$id,leg$id37)
    xx <- xx[!is.na(m),]
    w <- which(leg$id37 %in% xx$id)
    cat(w,file="lines.txt",sep="\n")
    
    message(hapfile,": reading ",length(w)," lines")
    
    haps <- fread(cmd=paste0("zcat ",hapfile," |",
                             "nl -nln -w1 -s' ' |",
                             "grep -F -w -f ~/scratch/gtex/lines.txt |",
                             "cut -d' ' -f ",cuts.eur),
                  fill=TRUE,sep=" ")
    haps  %<>% as.matrix()  %>% t()
    haps <- matrix(as(haps,"numeric"),nrow(haps),ncol(haps))
    colnames(haps) <- leg[w,]$id37

    vars <- apply(haps,2,var)
    drop <- which(vars==0)
    if(any(drop)) {
        haps <- haps[,-drop]
        xx <- xx[id %in% colnames(haps),]
    }
    
    ## finemap
    message("counting finemaps")
    xs <- split(xx,xx$gene_id)
    cres <- lapply(xs, function(xj) {
        hj <- haps[,xj$id,drop=FALSE]
        LD <- cor(hj,use="pair") ## std dev = 0 - need to drop 0 var variants at start
        if(nrow(xj)==1) { structure(1,names=xj$id) } else {
                with(xj,#[abs(tss_distance)<w,],
                     coloc::finemap.indep.signals(list(beta=slope,
                                                       varbeta=slope_se^2,
                                                       N=369,
                                                       MAF=maf,
                                                       type="quant",
                                                       snp=id),
                                                  LD=LD,
                                                  maxhits=6,
                                                  pthr=1,
                                                  method="mask"))
        }
    })
    names(cres) <- names(xs)
    saveRDS(cres,file=outfile)
}


## })  %>% do.call("rbind",.)


##     x1 <- x[abs(tss_distance)<100000,]
##     x1 <- x1[order(gene_id,pval_nominal),]
##     x1[,i:=1:.N,by="gene_id"]
##     x1 <- x1[i==1,]
##     sign(x1$slope)  %>% table()

##     x[,minp:=min(pval_nominal)]
## lapply(xs, function(z) sign(z[pval_nominal==min(pval_nominal),]$slope)) %>%unlist()%>%table()
## .

################################################################################

setwd("~/scratch/gtex")
library(magrittr)
library(data.table)

## number of genes
ng <- readRDS("ng.rds")
if(FALSE) {
    system("wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz")
    ng <- fread(cmd="zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz")
    ng <- nrow(ng)
    saveRDS(ng,"ng.rds")
}

## calc average number of snps in each distance window
nmaf <- readRDS("nmaf.rds")
if(FALSE) {
    RET <- vector("list",length(windows))
    j <- 1
    for(ichr in 1:22) {
        legfile <- file.path(dref,paste0("1000GP_Phase3_chr",ichr,".legend.gz"))
        message("reading legend ",basename(legfile))
        leg <- fread(cmd=paste("zcat",legfile))
        leg <- leg[EUR>0.005,] # limit to MAF>0.005 b/c don't see any lead snps with smaller MAF than this
        leg[,d:=cut(position,seq(min(position),max(position)+windows[1],by=windows[1]),include.lowest=TRUE)] 
        leg[,maf:=cut(pmin(EUR,1-EUR),seq(0,0.5,by=0.05),include.lowest=TRUE)] 
        RET[[j]] <- leg[,(n=.N),by=c("maf")]
        ## c(w,tot1=sum(nw1$V1),n1=nrow(nw1),
        ##                   tot5=sum(nw5$V1),n5=nrow(nw5))
            j <- j+1
    }
    ret <- do.call("rbind",RET)  %>% as.data.table()
    nmaf <- ret[,.(N=sum(V1)),by="maf"]
    nmaf <- nmaf[order(maf),]
    saveRDS(nmaf,file="nmaf.rds")
}

## 373 / 100kb MAF=0.01
## 274 / 100kb MAF=0.05

## number of signif lead snps in successive windows
x=fread(cmd="zcat GTEx_Analysis_v7_eQTL/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz")
x[,c("chr","pos.37","a1","a2","b37"):=tstrsplit(variant_id,"_")]
x[,pos.37:=as.numeric(pos.37)]

files <- list.files(".",pattern="^map3.*")
files
## files <- sub("3","2",files)
m <- lapply(files,readRDS)  %>% do.call("c",.)
head(m)
n <- sapply(m,length)
m <- m[n>0]
n <- sapply(m,length)
dt <- data.table(gene_id=rep(names(m),times=n),
                 id=unlist(lapply(m,names)),
                 z=unname(unlist(m)))
x[,id:=sub("_b37","",variant_id)]
x[,id:=gsub("_",":",id)]

head(dt)
head(x)
dt <- merge(dt,x[,.(id,gene_id,tss_distance,maf)],by=c("id","gene_id"))
dt[,cmaf:=cut(maf,breaks=seq(0,0.5,by=0.05),include.lowest=TRUE)]

par(mfrow=c(1,2))
hist(dt$maf,breaks=10,main="MAF of lead SNPs in GTeX")
barplot(nmaf$N,names.arg=levels(nmaf$maf),main="MAF of all SNPs in 1000 Genomes, EUR")
dt[,tss:=abs(tss_distance)]

## final - sigsnps/totalsnps
RET <- vector("list",length(windows))
for(i in seq_along(windows)) {
    w <- windows[i]
    tmp <- dt[tss<w,.(w=w,nsig=.N),by="cmaf"]
    tmp <- merge(tmp,nmaf,by.x="cmaf",by.y="maf")
    tmp[,N:=N*i] # adjust for window size
    tmp[,psig:=nsig/N]
    RET[[i]] <- tmp
}
ret <- rbindlist(RET)
ret[,newmaf:=(as.numeric(cmaf)+1) %/% 2]
ret[,newmaf:=factor(newmaf)]
levels(ret$newmaf) <- c("[0,0.1)","-0.2","-0.3","-0.4","(0.4,0.5]")
ret <- ret[,.(nsig=sum(nsig),N=sum(N)),by=c("newmaf","w")]
ret[,psig:=nsig/N]

library(ggplot2)
library(cowplot)
library(ggpubr)
library(viridisLite)
library(ggrepel)
## ggplot(ret,aes(x=cmaf,y=psig,group=factor(w),col=factor(w))) + geom_point() + geom_path()

cols=c("firebrick","DarkCyan")
p <- ggplot(ret,aes(x=w/1000,y=psig,group=newmaf,col=newmaf)) +
  geom_path() +
  geom_point(col="black") +
  geom_label(aes(label=newmaf),data=ret[w==windows[1] & as.numeric(newmaf) %in% c(1,5),],size=3,hjust=0.7) +
  ## annotate(x=windows[1]/1000,y=1.2e-4,label="q==10^{-4}",geom="label",parse=TRUE, col="darkblue",size=2.5) +
  geom_hline(yintercept=1e-4,linetype="dashed") +
  labs(x="Distance from TSS (kb)",
       y=expression(q~x~10^{4})) +
  theme_pubr() +
  theme(#axis.text=element_text(size=8),
        legend.position="none") +
  ## scale_colour_manual(values = rainbow(10)) +
  ## scale_color_gradient(low = "blue", high = "red") +
  scale_y_continuous(breaks=c(0.0001,0.0005,0.001),labels=c(1,5,10),limits=c(5e-5,0.001))+
  ## scale_y_log10(breaks=c(0.0001,0.0003,0.001),labels=c(1,3,10)) +
  xlim(-50,1000) +
  scale_colour_viridis_d("MAF",begin=0,end=0.5, guide = guide_legend(reverse = TRUE)) +
  background_grid()
plot_grid(p,labels=c("a",""),ncol=1)

w <- 1.2
ggsave("~/fig-coloc-gtex.pdf",height=w*4,width=w*8/3)

