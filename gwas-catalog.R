library(magrittr)
library(data.table)
## x=fread("curl -s https://www.ebi.ac.uk/gwas/api/search/downloads/alternative")
## save(x,file="~/gwascat.RData")
(load("~/gwascat.RData"))

setnames(x, make.names(names(x)))
    x <- x[GENOTYPING.TECHNOLOGY=="Genome-wide genotyping array",]
    x[,ccqtl:=ifelse(grepl("cases",INITIAL.SAMPLE.SIZE),"Case/control","Quantitative trait")]
x[,INITIAL.SAMPLE.SIZE:=gsub(", ",", ",INITIAL.SAMPLE.SIZE)]
x[,INITIAL.SAMPLE.SIZE:=gsub(" and ",", ",INITIAL.SAMPLE.SIZE)]
## x[grep("protein",DISEASE.TRAIT),ccqtl:="pqtl"] 
    ss <- strsplit(x$INITIAL.SAMPLE.SIZE,", ")  %>% lapply(., function(x) gsub(" .*|,","",x))  %>% lapply(.,as.numeric)  %>% lapply(.,sum)  %>% unlist()
    head(ss)
    x[,ntot:=ss]

    ss <- strsplit(x$INITIAL.SAMPLE.SIZE,", ")  %>% lapply(., function(x) x[grep("cases",x)])  %>% lapply(., function(x) gsub(" .*|,","",x))  %>% lapply(.,as.numeric)  %>% lapply(.,sum)  %>% unlist()
    head(ss)
    x[,ncases:=ss]

    ss <- strsplit(x$INITIAL.SAMPLE.SIZE,", ")  %>% lapply(., function(x) x[grep("controls",x)])  %>% lapply(., function(x) gsub(" .*|,","",x))  %>% lapply(.,as.numeric)  %>% lapply(.,sum)  %>% unlist()
    head(ss)
    x[,ncontrols:=ss]
    x[ncontrols==0 & ccqtl=="Case-control",ccqtl:="Quantitative Trait"]

    head(x)

    ## remove outliers
    x <- x[!(STUDY %in% c("The coexistence of copy number variations (CNVs) and single nucleotide polymorphisms (SNPs) at a locus can result in distorted calculations of the significance in associating SNPs to disease.")),]
with(x[ccqtl=="Case/control" & P.VALUE<5e-8], hist(abs(log(OR.or.BETA)),xlim=c(0,10),breaks=100))

    
    nx <- x[P.VALUE<5e-8, .(hits=.N), by=c("STUDY","DISEASE.TRAIT","MAPPED_TRAIT","P.VALUE..TEXT.","ccqtl","ntot","ncases","ncontrols")]
    nx[ccqtl=="Case/control",ntot:=ncases]
    head(nx)
    ## classify traits



tt <- table(nx$DISEASE.TRAIT)
    table(tt)
    tt[tt>100]
    tt <- table(nx$STUDY)
    table(tt)
    tt[tt>100]

    library(ggplot2)
    library(cowplot)
    with(nx,cor.test(ntot,hits,use="pair"))

ann <- data.frame(y=c(100,100,1000),x=c(1e+5,1e+5,1e+6),ccqtl=c("Case/control","Quantitative trait","Quantitative trait"),
                  label=c("10,000 samples\n100 hits",
                          "10,000 samples\n100 hits",
                          "100,000 samples\n1000 hits?"))
library(ggrepel)
library(ggpubr)
nx[,n:=ifelse(ccqtl=="Quantitative trait",ntot,pmin(ncases,ncontrols))]
nx <- nx[!is.na(n),]
nx[,hmin:=min(hits,na.rm=TRUE),by="ccqtl"]
nx[,hmax:=max(hits,na.rm=TRUE),by="ccqtl"]
nx <- nx[order(n,decreasing=TRUE),]
ann <- nx[hits==hmin | hits==hmax,]
ann <- unique(ann,by=c("hits","ccqtl"))
ann

p=ggplot(nx[n>0,],aes(x=(n),y=(hits))) + geom_point() +
                                        ## geom_smooth(method="lm") +
  facet_wrap(~ccqtl) + scale_x_log10() + scale_y_log10() + labs(x="Total samples (log scale)",y="Number of hits (log scale)") +
    geom_label_repel(aes(label=DISEASE.TRAIT),data=ann,hjust=1.2,vjust=-1,#fontface="bold",
                     col="darkblue",size=2.5) +
  ## geom_hline(aes(yintercept=y),data=ann,linetype="dashed",col="darkorange") +
  ## geom_vline(aes(xintercept=x),data=ann,linetype="dashed",col="darkorange") +
  ## geom_label(aes(x=x,y=y,label=label),data=ann,hjust=1,fill="darkorange",colour="white",fontface="bold") +
    theme_pubr()
library(cowplot)
plot_grid(p,labels="b")
ggsave("gwas-catalog.pdf",height=4,width=2*8/3)

ann


##   #annotate("label",label="10,000 samples\n100 hits", x = 1e+5, y = 1e+2,
##                                                                                                                                       vjust = 0, hjust = 0, fontface = 'bold')
## #geom_label(x=1e+5,y=100,label="10,000 samples\n100 hits") #c(100,500))

## ggplot(nx,aes(x=(ntot),y=(hits))) + geom_point() +
##                                         #geom_smooth() +
##   facet_wrap(~ccqtl) + scale_x_log10() + scale_y_log10() + labs(x="Total samples (log scale)",y="Number of hits (log scale)") +
##   geom_hline(yintercept=1000,linetype="dashed",col="orange") + geom_vline(xintercept=1e+6,linetype="dashed",col="orange") +
##   annotate("label",label="100,000 samples\n1000 hits?", x = 1e+6, y = 1e+3,
##                                                                                                                                       vjust = 1, hjust = 1, fontface = 'bold')
## #geom_label(x=1e+5,y=100,label="10,000 samples\n100 hits") #c(100,500))

##     nx[hits>900,]


##     ## f <- function(beta,se,pr) {
##     ##     w <- 0.2
##     ##     lABF
##     library(coloc)
##     x[,RISK.ALLELE.FREQUENCY:=as.numeric(RISK.ALLELE.FREQUENCY)]
## x <- x[RISK.ALLELE.FREQUENCY<1,] # because - MAF=3?!?
## x[,MAF:=pmin(RISK.ALLELE.FREQUENCY,1-RISK.ALLELE.FREQUENCY)] # because - MAF=3?!?
## x[,P.VALUE:=as.numeric(P.VALUE)]
##     x[ccqtl=="Quantitative trait",bf:=coloc:::approx.bf.p(p=P.VALUE,f=RISK.ALLELE.FREQUENCY,type="quant",N=ntot)[,"lABF"]]
##     x[ccqtl=="Case-control",bf:=coloc:::approx.bf.p(p=P.VALUE,f=RISK.ALLELE.FREQUENCY,type="cc",N=ntot,s=ncases/ntot)[,"lABF"]]

## ggplot(x,aes(x=bf,y=-log10(P.VALUE))) + geom_point()

## pr <- 1e-4
## x[,pp4:=pr * exp(bf)/(1-pr+pr*exp(bf))]
## pr <- 5e-4
## x[,pp54:=pr * exp(bf)/(1-pr+pr*exp(bf))]
## pr <- 5e-5
## x[,pp55:=pr * exp(bf)/(1-pr+pr*exp(bf))]
## m <- melt(x[,.(P.VALUE,pp4,pp54,pp55,ntot,MAF,ccqtl)],c("ccqtl","MAF","ntot","P.VALUE"))

## library(ggpubr)
## xx <- x[!is.na(MAF)&ntot>1000 & !is.na(PVALUE_MLOG) & !is.na(pp4) & MAF!=0,]
## ggplot(xx,aes(x=PVALUE_MLOG,y=pp4,group=cut(MAF,c(0,0.05,0.5)))) +
##                                          ## log10(ntot))) +
##   geom_point(aes(col=cut(MAF,c(0,0.05,0.5))),size=1,alpha=0.3) + 
##   geom_smooth(col="darkolivegreen",se=FALSE,data=xx[MAF>0.05,]) +
##   facet_wrap(~ccqtl) +
##   scale_colour_manual("MAF",values=c("grey30","skyblue")) +#scale_colour_gradient(low="#cccccc",high="#000000")+
##   geom_vline(xintercept=-log10(5e-8),linetype="dashed") + geom_hline(yintercept=0.9,linetype="dashed") +
##   scale_x_continuous(breaks=seq(6,14,2),limits=c(5,15),
##                      labels=sapply(seq(6,14,2),function(z) bquote(10^.(z)))) +
##   scale_y_continuous(breaks=seq(0,1,0.2)) +
##   background_grid() + labs(x="-log10 p value",y="Posterior probability") +
##   theme_pubr() +
##   theme(legend.position=c(0.9,0.1),legend.background = element_rect(colour="black"))


## xx[,cMAF:=cut(MAF,c(0,0.05,0.5))]
## ggscatter(xx, x = "PVALUE_MLOG", y = "pp4",
##                 add = "loess",               # Add regression line
##                 conf.int = FALSE,                # Add confidence interval
##           facet.by="ccqtl",
##           color = "cMAF",
##           palette = "jco", # Color by groups "cyl"
##                 shape = "cMAF"                   # Change point shape by groups "cyl"
##                 )+
##   stat_cor(aes(color = cMAF), label.x = 3)       # Add correlation coefficient
