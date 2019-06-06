#!/usr/bin/env Rscript
library(randomFunctions)
library(snpStats)
library(magrittr)
devtools::load_all("~/RP/coloc")
source("~/DIRS.txt")
args <- getArgs(defaults=list(N=1000,NSIM=100,NCV=3,SPECIAL=2),
                numeric=c("N","NSIM","NCV"))

d <- COLOCINDEP

## read in haplotypes
library(data.table)
patt <- "cvaryr2-v13"
files <- list.files(d,pattern=patt,full=TRUE)
message("files found: ",length(files))

data <- lapply(files, function(f) {
    ## cat(f,"\t")
    x=eval(as.symbol(load(f)))
    cols <- c(grep("PP",names(x),value=TRUE),"nsnps")
    x[ ,(cols) := lapply(.SD, unlist), .SDcols = cols] 
    x[,f:=basename(f)]
    x})  %>% rbindlist(.,fill=TRUE)
data[,group:=factor(paste(f,isim))]
table(data$method)
table(data$NCV,data$special)
data <- unique(data)

library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggridges)

source("~/Projects/coloc-cond-mask/plotspec.R")

## pdf("~/kk.pdf",width=8,height=6)
## data[,H4.H3:=PP.H4.abf/PP.H3.abf]

## significance

## assign "hits"
hit.r2thr <- 0.8
results <- copy(data)
## results <- results[r2.finemap %in% c(NA,0.01),]
results[,minz:=min(pmin(abs(hit1.margz),abs(hit2.margz))),by=c("f","isim","method","r2.finemap")]
results[,pwr:=ifelse(minz < 9.35,"low","high")]
## results <- results[!is.na(r2.cvA1.hits),] # previous bad script included cond dt twice
results[,tested.cv1:="?"][,tested.cv2:="?"]
results[r2.cvA1.hits > r2.cvB1.hits & r2.cvA1.hits > hit.r2thr, tested.cv1:="A"]
results[r2.cvB1.hits > r2.cvA1.hits & r2.cvB1.hits > hit.r2thr, tested.cv1:="B"]
results[r2.cvA2.hits > r2.cvB2.hits & r2.cvA2.hits > hit.r2thr, tested.cv2:="A"]
results[r2.cvB2.hits > r2.cvA2.hits & r2.cvB2.hits > hit.r2thr, tested.cv2:="B"]
## results[pmax(r2.cvA1.hits,r2.cvB1.hits)>hit.r2thr,tested.cv1:=ifelse(r2.cvA1.hits > r2.cvB1.hits, "A", "B")]
results[NCV=="onetwo" ,tested.cv2:="A"]
results[NCV=="onetwo" & special==0,tested.cv2:="C"]
## results[NCV=="twotwo" & special==0 & tested.cv2=="A", tested.cv2:="C"]
## remap A, B -> C, D for trait two, then map back
results[NCV=="twotwo" & tested.cv2=="B", tested.cv2:="D"]
results[NCV=="twotwo" & tested.cv2=="A", tested.cv2:="C"]
## remap C->A for 1-5
results[NCV=="twotwo" & special %in% c(1,2,3,4,5) & tested.cv2=="C",tested.cv2:="A"]
## remap D->B for 2,4
results[NCV=="twotwo" & special %in% c(2,4) & tested.cv2=="D",tested.cv2:="B"]
## remap D->C for 1,3,5
results[NCV=="twotwo" & special %in% c(1,3,5) & tested.cv2=="D",tested.cv2:="C"]
with(results,table(paste(tested.cv1,tested.cv2,sep="-"),method,group=paste(NCV,special,sep="-")))
results[,tested.cv:=paste(tested.cv1,tested.cv2,sep=":")]
setnames(results,
         grep("PP",names(results),value=TRUE),
         gsub("PP.|.abf","",grep("PP",names(results),value=TRUE)))
results[,special:=paste(NCV,special,sep="-")]

## options(digits=2)
## results[special=="twotwo-4" & f=="cvaryr2-v10243847a51be35.RData" & r2.finemap %in% c(NA,0.01) & method %in% c("single","cond","mask") & isim==2,.(method,tested.cv,hit1,hit2,H3,H4, hit1.margz ,hit2.margz)]
## results[special=="twotwo-4" & f=="cvaryr2-v10243847a51be35.RData" & isim==3 & r2.finemap %in% c(NA,0.01) & method %in% c("single","cond","mask") ,.(method,tested.cv,hit1,hit2,H3,H4, hit1.margz ,hit2.margz,r2.cvA1.hits,r2.cvA2.hits,r2.cvB1.hits,r2.cvB2.hits,r2.tr1)]

## sresults <- results[,.(h0=sum(H0),
##                        h1=sum(H1) + sum(H2),
##                        h3=sum(H3),
##                        h4=sum(H4)),by=c("f","isim","method","special")]
## m <- melt(sresults[,.(special,method,h0,h1,h3,h4)],
##           c("special","method"))
## m[,pc:=value/sum(value),by=c("special","method")]
  library(viridis)
library(ggplot2)
library(ggpubr)
## ggplot(m,aes(x=method,y=pc,fill=variable)) + geom_col() +
##   facet_grid(special ~.) + theme_pubr() +
## scale_fill_viridis_d()
## next: look at best1, best2 and the CVs in details.
## 1. are the bests in LD with the CVs?
## 2. are

## simulations options, NCV=onetwo
## 2 = share 1, weaker effect
## 1 = share 1, equal effect (strongest), + indep each
## 0 = share 0
## indep want H4=1, H3=1, H1=1, H2=1
## step want H4=1, H3=3?
## 0 = share 0

## total conclusions
library(cowplot)
byvars <- c("special","method","tested.cv","r2.finemap","r2.process") #,"ld")
byvars <- setdiff(byvars,"tested.cv")
kk <- results[abs(hit1.margz)>4.89 & abs(hit2.margz)>4.89#NCV=="onetwo"
            , .(h0=sum(H0),h1=sum(H1),h2=sum(H2),h3=sum(H3),h4=sum(H4)),by=c(byvars,"f","isim")]
kk <- kk[, .(ndata=length(unique(paste(f,isim))),N=.N,h0=sum(h0),h1=sum(h1),h2=sum(h2),h3=sum(h3),h4=sum(h4)), by=c(byvars)]
kk[,c("h0","h1","h2","h3","h4"):=list(h0/ndata,h1/ndata,h2/ndata,h3/ndata,h4/ndata)]
m <- melt(kk,c(byvars,"N","ndata"))
## m <- m[method %in% c("single","cond","mask","condmask","maskcond") & (is.na(r2.finemap) |r2.finemap==0.01),]
m <- m[method %in% c("single","cond","mask") , ]
## m[,method:=factor(method,levels=c("single","cond","condmask","mask"))]
m[,x:="x"] #paste(ld)]
p.total <- ggplot(m,#!grepl("?",tested.cv,fixed=TRUE) &
## (r2.finemap==0.01|method %in% c("s),],
                  aes(x=x,#paste(method,r2.finemap),
                      y=value,fill=variable)) + geom_col() +
  facet_grid(special ~ method + r2.finemap ,space="free_x",scales="free_x") + theme_pubr() +
  ## facet_grid(special ~ ld,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hypoth.") +
  ylab("Avg. posterior") +
  xlab("Method") +
  background_grid() +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        ## strip.text.y = element_blank(),
        axis.text.x = element_text(angle=-90)
        )
p.total

m[special=="twotwo-4",]
m[special=="twotwo-4",sum(value),by=byvars]

## number tested - need r2.finemap close to 0.1 to get similar numbers to conditional
byvars <- c("special","method","tested.cv","r2.finemap","r2.process","ld")
results[,ndata:= length(unique(paste(f,isim))), by=setdiff(byvars,"tested.cv")]
kk <- results[abs(hit1.margz)>4.89 & abs(hit2.margz)>4.89#NCV=="onetwo"
            , .(N=.N,h0=sum(H0),h1=sum(H1),h2=sum(H2),h3=sum(H3),h4=sum(H4)), by=c("ndata",byvars)]
kk[,c("h0","h1","h2","h3","h4"):=list(h0/ndata,h1/ndata,h2/ndata,h3/ndata,h4/ndata)]
## , .(N=.N,h0=sum(H0)/ndata,h1=sum(H1)/ndata,h2=sum(H2)/ndata,h3=sum(H3)/ndata,h4=sum(H4)/ndata), by=byvars]
## , .(N=.N,h0=mean(H0),h1=mean(H1),h2=mean(H2),h3=mean(H3),h4=mean(H4)), by=byvars]
## kk <- results[NCV=="onetwo", .(H0=sum(H0),H1=sum(H1),H2=sum(H2),H3=sum(H3),H4=sum(H4)), by=c("special","method","tested.cv")]
## kk[method=="cond",r2.finemap:=0]
## p.r <- ggplot(kk[#!grepl("?",tested.cv,fixed=TRUE) &
## (r2.process==0.01|method=="cond"),],
##               aes(x=tested.cv,y=N,fill=factor(r2.finemap))) + geom_col(position="dodge") +
##   facet_grid(special ~ tested.cv,space="free_x",scales="free") + theme_pubr() +
##   ## scale_fill_viridis_d("hypoth.") +
##   ylab("Avg. posterior") +
##   xlab("Tested variants") +
##   theme(legend.position="right",
##         panel.border=element_rect(linetype="solid",fill=NA),
##         ## strip.text.y = element_blank(),
##         axis.text.x = element_text(angle=-90)
##         )
## p.r

## vary r2
m <- melt(kk,c(byvars,"N","ndata"))
m[,method:=factor(method,levels=unique(data$method))] #c("single","cond","condmask","mask"))]
m <- m[method %in% c("single","cond","mask","condmask","maskcond") ## & ld=="lowld"
       & !grepl("?",tested.cv,fixed=TRUE),]
## p.r <- ggplot(m[is.na(r2.finemap) | r2.finemap %in% c(0.01,0.02,0.1),],#!grepl("?",tested.cv,fixed=TRUE) &
p.r <- ggplot(m[method %in% c("single","cond","mask"),],#!grepl("?",tested.cv,fixed=TRUE) &
## (r2.finemap==0.01|method %in% c("s),],
              aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special  ~ method + r2.finemap ,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hyp.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid() +
  ylim(0,1.05) + 
  theme_pubr() +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        ## strip.text.y = element_blank(),
        axis.text.x = element_text(angle=-90)
        )
p.r

## pub plots
m <- melt(kk,c(byvars,"N","ndata"))
m[grepl("?",tested.cv,fixed=TRUE),tested.cv:="?"]
m[,method:=factor(method,levels=unique(data$method))] #c("single","cond","condmask","mask"))]
m[,ncv:=sub("-.*","",special)]
m[,r2:=ifelse(method=="cond","cond",r2.finemap)]

L <- lapply(0:3, function(spec) {
 ggplot(m[method %in% c("cond","mask","single") & special==paste0("onetwo-",spec) & r2.finemap %in% c(NA,0.01),], aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special  ~ method ,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hyp.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid() +
  scale_y_continuous(breaks=c(0,1),limits=c(0,1)) +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=-90)
        )
})
p.1 <- plot_grid(plotlist=L,ncol=1)
pp <- plot_grid(spec.1,p.1,nrow=1,rel_widths=c(0.2,0.8))
w <- 1.2
save_plot("~/fig-coloc-onetwo.pdf",pp,base_height=w*8,base_width=w*6)


L <- lapply(0:3, function(spec) {
 ggplot(m[special==paste0("onetwo-",spec) & r2.finemap %in% c(NA,0.01),], aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special  ~ method ,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hyp.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid() +
  scale_y_continuous(breaks=c(0,1),limits=c(0,1)) +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=-90)
        )
})
p.1 <- plot_grid(plotlist=L,ncol=1)
pp <- plot_grid(spec.1,p.1,nrow=1,rel_widths=c(0.2,0.8))
w <- 1.2
save_plot("~/fig-coloc-onetwo-mixed.pdf",pp,base_height=w*8,base_width=w*6)

L <- lapply(0:3, function(spec) {
 ggplot(m[special==paste0("onetwo-",spec) & method %in% c("cond","mask"),], aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special  ~ r2 ,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hyp.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid() +
  scale_y_continuous(breaks=c(0,1),limits=c(0,1)) +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=-90)
        )
})
p.1 <- plot_grid(plotlist=L,ncol=1)
pp <- plot_grid(spec.1,p.1,nrow=1,rel_widths=c(0.2,0.8))
w <- 1.2
save_plot("~/fig-coloc-onetwo-r2.pdf",pp,base_height=w*8,base_width=w*6)
## w <- 1.2
## save_plot("~/coloc-gsk.pdf",p,base_height=w*6,base_width=w*8)


L <- lapply(0:5, function(spec) {
 ggplot(m[method %in% c("cond","mask","single") & special==paste0("twotwo-",spec) & r2.finemap %in% c(NA,0.01),], aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special  ~ method ,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hyp.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid() +
  scale_y_continuous(breaks=c(0,1),limits=c(0,1)) +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=-90)
        )
 })
p.2 <- plot_grid(plotlist=L,ncol=1)
pp <- plot_grid(spec.2,p.2,nrow=1,rel_widths=c(0.2,0.8))
w <- 1.2
save_plot("~/fig-coloc-twotwo.pdf",pp,base_height=w*8,base_width=w*6)


L <- lapply(0:5, function(spec) {
 ggplot(m[special==paste0("twotwo-",spec) & r2.finemap %in% c(NA,0.01),], aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special  ~ method ,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hyp.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid() +
  scale_y_continuous(breaks=c(0,1),limits=c(0,1)) +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=-90)
        )
 })
p.2 <- plot_grid(plotlist=L,ncol=1)
pp <- plot_grid(spec.2,p.2,nrow=1,rel_widths=c(0.2,0.8))
w <- 1.2
save_plot("~/fig-coloc-twotwo-mixed.pdf",pp,base_height=w*8,base_width=w*6)


L <- lapply(0:5, function(spec) {
 ggplot(m[special==paste0("twotwo-",spec) & method %in% c("cond","mask"),], aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
  facet_grid(special  ~ r2 ,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hyp.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid() +
  scale_y_continuous(breaks=c(0,1),limit=c(0,1)) +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=-90)
        )
 })
p.2 <- plot_grid(plotlist=L,ncol=1)
pp <- plot_grid(spec.2,p.2,nrow=1,rel_widths=c(0.2,0.8))
w <- 1.2
save_plot("~/fig-coloc-twotwo-r2.pdf",pp,base_height=w*8,base_width=w*6)


















## simplified version for talk
kk <- results[method!="cond.manual" & !is.na(H0) & NCV=="onetwo" & tested.cv1!="?", .(h0=mean(H0),h1=mean(H1),h2=mean(H2),h3=mean(H3),h4=mean(H4)), by=c("special","method","tested.cv")]
## kk <- results[!is.na(H0) & NCV=="onetwo", .(h0=sum(H0),h1=sum(H1),h2=sum(H2),h3=sum(H3),h4=sum(H4)), by=c("special","method","tested.cv")]
m <- melt(kk,c("special","method","tested.cv"))
m[,method:=factor(method,levels=c("single","cond.manual","cond","mask"))]
L <- lapply(0:2, function(spec) {
    ggplot(m[special==spec,],aes(x=tested.cv,y=value,fill=variable)) + geom_col() +
      facet_grid(. ~ method,space="free",scales="free") + theme_pubr() +
      scale_fill_viridis_d("hypoth.") +
      ylab("Avg. posterior") +
      xlab("Tested variants") +
      theme(legend.position="right",
            panel.border=element_rect(linetype="solid",fill=NA),
            strip.text.y = element_blank()
            )
})
p.r <- plot_grid(plotlist=L,ncol=1)



plot_grid(p.l,p.r,nrow=1,rel_widths=c(0.2,0.8))



pp <- plot_grid(p.l,p.r,nrow=1,rel_widths=c(0.2,0.8))

w <- 1.2
## save_plot("~/coloc-gsk.pdf",p,base_height=w*6,base_width=w*8)

## results[,tested.cv1:=0]
data <- copy(results)

## how often do hits change, with r2?
y <- data[,.(nhits1=length(unique(hit1)),
    nhits2=length(unique(hit2))),
    by=c("group","method","NCV","ntested","special")]
head(y)
yy <- y[,.(count=.N),by=c("method","NCV","special","nhits1","nhits2")]
yy <- yy[order(count,decreasing=TRUE),]
head(yy)

## ggplot(y,aes(x=special,y=nhits1)) + geom_violin() + facet_wrap(method ~ NCV)

with(y[method=="mask",],table(nhits1,nhits2,NCV,special))
with(y[method=="mask",],table(ntested,NCV,special))


## number of things tested
## increases with r2.finemap
## largely indep of r2.process - as expected
## ? indep of ld structure
library(viridis)
cols <- grep("PP",names(data))
m <- data[,.(value=mean(ntested,na.rm=TRUE)), by=c("NCV","special","method","r2.finemap","r2.process","ld")]
head(m)
m[,special:=paste0("CV-",NCV,"-",special)]
setnames(m,"special","sim")

table(m$r2.process)
m[is.na(r2.finemap),r2.finemap:=0]
mc <- m[method=="cond",]
m <- m[method!="cond",]
setnames(mc,"value","yint")
m <- merge(m,mc[,.(sim,ld,yint)],by=c("sim","ld"))

ggplot(m, aes(x=r2.finemap,y=value,col=factor(r2.process))) +
  geom_point() +#range() +
  geom_path() +
  geom_hline(aes(yintercept=yint)) +
  ## geom_col() +
  ## scale_x_log10(breaks=c(1e-7,1e-5),labels=c("-7","-5")) +
  background_grid() +
  scale_y_continuous(limits=c(0,8)) +
  ## scale_y_log10() +
  ## scale_colour_viridis_d("post. hyp") +
  ## scale_fill_viridis_d("post. hyp") +
  ## facet_grid(sim~ r2.process,labeller = label_both) +
  facet_grid(sim~ ld,labeller = label_both) +
  theme_pubr() +
  labs(x="r2.finemap",y="ntested") +
  theme(legend.position="right")


## average conclusion
## r2.process needs to be small - 0.01
library(viridis)
cols <- grep("PP",names(data))
m <- melt(data,measure.vars=cols)
m <- m[,.(value=mean(value)), by=c("NCV","special","method","ntested","r2.finemap","r2.process","variable")]
head(m)
m[,special:=paste0("CV-",NCV,"-",special)]
m[,variable:=gsub("PP.|.abf","",variable)]
setnames(m,"special","sim")
m[is.na(r2.finemap),r2.finemap:=0]
m[is.na(r2.process),r2.process:=0]

table(m$ntested)
table(m$r2.process)
m <- m[ntested %in% c(1,2,4),]

ggplot(m[(method=="cond" | r2.finemap==0.1),], aes(x=r2.process,y=value,fill=variable,col=variable)) +
  ## geom_violin() +#range() +
  geom_col() +
  ## scale_x_log10(breaks=c(1e-7,1e-5),labels=c("-7","-5")) +
  background_grid() +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  ## scale_y_log10() +
  scale_colour_viridis_d("post. hyp") +
  scale_fill_viridis_d("post. hyp") +
  facet_grid(sim~ ntested,labeller = label_both) +
  theme_pubr() +
  labs(x="r2.process",y="mean posterior") +
  theme(legend.position="right")




m <- m[p12<1  & special %in% c("3","4"), ]
m[,group:=paste0(N,special)]
## m <- rbind(m[special=="3",],m[special=="4",])
eps <- 1e-16
m[H4.H3<eps,H4.H3:=eps]
m[H4.H3>1/eps,H4.H3:=1/eps]
m <- melt(m[,.(N,p12,PP.H3.abf,PP.H4.abf,special)],c("N","p12","special"))
m[,special:=paste0("H",special)]
m <- m[order(p12),]

ggplot(m, aes(x=value,y=special,fill=variable,col=variable)) +
  geom_density_ridges(alpha=0.5) +
  facet_grid(p12 ~ N)

pdf("~/coloc-varyprior-v2.pdf",height=8,width=8*4/3)
library(viridis)
m <- melt(data[p12<1,],c("p12","N","special"),grep("abf",names(data)))
m <- m[,.(value=mean(value)),by=c("N","p12","special","variable")]
m[,special:=paste0("H",special)]
m[,variable:=gsub("PP.|.abf","",variable)]
setnames(m,"special","sim")
ggplot(m, aes(x=N,y=value,fill=variable,col=variable)) +
  ## geom_violin() +#range() +
  geom_col() +
  scale_x_log10(breaks=c(500,5000)) +
  ## scale_x_log10(breaks=c(1e-7,1e-5),labels=c("-7","-5")) +
  background_grid() +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  ## scale_y_log10() +
  scale_colour_viridis_d("post. hyp") +
  scale_fill_viridis_d("post. hyp") +
  facet_grid(sim~ p12,labeller = label_both) +
  theme_pubr() +
  labs(x="Sample size",y="mean posterior") +
  theme(legend.position="right")
dev.off()


## pdf("~/coloc-varyprior.pdf",height=6,width=8)
m <- data[p12<1,.(PP.H4.abf=mean(PP.H4.abf,na.rm=TRUE),l=quantile(PP.H4.abf,0.05),u=quantile(PP.H4.abf,0.95)),by=c("p12","N","special")]
m[,special:=paste0("H",special)]
m <- m[order(N),]
setnames(m,"special","hyp")
ggplot(m, aes(x=N,y=PP.H4.abf,ymin=l,ymax=u,fill=hyp,col=hyp)) +
  ## geom_violin() +#range() +
  geom_path() +
  geom_point() +
  scale_x_log10(breaks=c(500,5000)) +
  background_grid() +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  ## scale_y_log10() +
  facet_grid(hyp~ p12,labeller = label_both) +
  theme_pubr() +
  labs(x="Sample size",y="median expected posterior, H4") +
  theme(legend.position="none")
## dev.off()

x <- melt(data[isim<=20,c(grep("PP",names(data),value=TRUE),"special","p12","N","group"),with=FALSE],
          c("special","p12","N","group"))
ggplot(x,aes(x=p12,y=value,group=group,col=factor(N))) + geom_path(alpha=0.1) + facet_wrap(special ~ variable) + scale_x_log10()


xm <- x[,.(value=mean(value)),by=c("special","p12","N","variable")]
ggplot(xm,aes(x=p12,y=value,group=factor(N),col=factor(N))) + geom_path() + facet_wrap(special ~ variable) + scale_x_log10()


       
x <- melt(data[,c(grep("PP",names(data),value=TRUE),"special","p12","N"),with=FALSE],
          c("special","p12","N"))
x[,variable:=gsub("PP.|.abf","",variable)]
x[,special:=as.integer(special)]
x[,p12:=as.factor(unlist(p12))]
x[,value:=unlist(value)]


## y <- x[,.(mv=mean(value)),by=c("special","p12","variable")]
y <- x[,.(mv=mean(value)),by=c("special","p12","variable","N")]
ggplot(y,aes(x=special,y=mv,fill=variable)) + geom_col() + facet_grid(N~p12)

library(ggridges)
ggplot(x, aes(x=value,y=factor(special),fill=variable,col=variable)) + geom_density_ridges() +
  facet_grid(N~p12)


## special=0..4 equiv to H_i

## probs <- function(nsnps,p12,p1=1e-4,p2=1e-4) {
##     h1 <- p1 * nsnps
##     h2 <- p2 * nsnps
##     h3 <- p1 * p2 * nsnps * (nsnps-1)
##     h4 <- p12 * nsnps
##     h0 <- 1 - h1 - h2 - h3 - h4
##     list(h0=h0,h1=h1,h2=h2,h3=h3,h4=h4)
## }
## data[,c("h0","h1","h2","h3","h4"):=probs(nsnps,p12)]

weights <- function(special,nsnps,p12,p1=1e-4,p2=1e-4) {
    h1 <- p1 * nsnps
    h2 <- p2 * nsnps
    h3 <- p1 * p2 * nsnps * (nsnps-1)
    h4 <- p12 * nsnps
    h0 <- 1 - h1 - h2 - h3 - h4
    h <- cbind(h0,h1,h2,h3,h4)
    h[ cbind(1:nrow(h) , as.numeric(special)+1) ] 
}
for(ip in unique(data$p12)[4]) {
    thing <- paste0("w.",ip)
    data[,(thing):=weights(special,nsnps,ip)]
}
data[,w:=NULL]
setnames(data,"w.5e-06","w")
## data[,w:=weights(N,special,nsnps,p12)]

## data[,right:=ifelse(special=="0",PP.H0.abf,
##              ifelse(special=="1",PP.H1.abf,
##              ifelse(special=="2",PP.H2.abf,
##              ifelse(special=="3",PP.H3.abf,
##                              PP.H4.abf))))]
## data[,right0:=right/(1-PP.H0.abf)]

library(ggplot2)
x <- melt(data[,c("special","N",
                  grep("PP.",names(data),value=TRUE),
                  grep("^w",names(data),value=TRUE),"p12"),with=FALSE],
          c("special","p12","N",
            grep("PP.",names(data),value=TRUE)),
          variable.name="weight.p12",
          value.name="weight")
ypr <- x[,list(pr=sum(weight/N)),by=c("special","weight.p12","p12")]
ypr <- unique(ypr[,.(special,weight.p12,pr)])
ypr[,sum(pr),by="weight.p12"] # should be all 1
x[,special:=NULL]
x <- melt(x,
          c("N","p12","weight.p12","weight"),
          variable.name="hyp",
          value.name="pp")
x[,special:=gsub("PP.H|.abf","",hyp)]
ypp <- x[,list(pp=sum(pp*weight/N)),by=c("weight.p12","p12","special")]
ypp[,sum(pp),by=c("weight.p12","p12")] # should be all 1

y <- merge(ypr,ypp,by=c("weight.p12","special"))
tmp=y[weight.p12=="w.1e-04" & p12 %in% c(1e-04,1e-08),]
tmp

ggplot(y[weight.p12=="w.1e-04",],aes(x=pr,y=pp)) + geom_abline() + geom_point(aes(colour=special)) + facet_wrap(~p12)

ggplot(y,aes(x=pr,y=pp)) + geom_abline() + geom_point(aes(colour=special)) + facet_grid(weight.p12~p12)

y <- x[special!="0",list(mv=sum(right*value)/sum(value)),by=c("variable","p12")]
ggplot(y,aes(x=variable,y=mv,col=factor(p12),group=factor(p12))) + geom_point() + geom_path() + facet_wrap(~p12)


y0 <- x[special!="0",list(mv=sum(right0*value)/sum(value)),by=c("variable","p12")]
ggplot(y0,aes(x=p12,y=mv,col=variable)) + geom_point() + geom_path()




y0 <- data[special!="0",list(mv=sum(right0*w)/sum(w)),
          by=c("p12")]

ggplot(y,aes(x=p12,y=mv)) + geom_point() + ylim(0.5,1)
ggplot(y0,aes(x=p12,y=mv)) + geom_point() + ylim(0.5,1)

y <- data[,.(mv=mean(right)),by=c("p12")]
ggplot(y,aes(x=p12,y=mv)) + geom_point() + ylim(0.5,1)

y <- x[,
x <- melt(data[,c("right","special","p12","w"),with=FALSE],
          c("special","p12","w"))
x[,variable:=gsub("PP.|.abf","",variable)]
x[,special:=as.integer(special)]
x[,p12:=as.factor(unlist(p12))]
x[,value:=unlist(value)]


## y <- x[,.(mv=mean(value)),by=c("special","p12","variable")]


data[,
y <- data[,.(mv=mean(PP.H4.abf)),by=c("p12")]
