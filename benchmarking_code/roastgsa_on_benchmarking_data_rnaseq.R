
####### load R packages
library(Biobase)
library(DESeq2)
library(roastgsa)
library(GSEABenchmarkeR)
library(EnrichmentBrowser)
library(preprocessCore)
library(hgu133plus2.db)


## R version 4.1.3 (2022-03-10)
## EnrichmentBrowser_2.24.2
## GSEABenchmarkeR_1.14.0
## preprocessCore_1.56.0
## GSEABenchmarkeR_1.2.1

####### create destination folder

####### create destination folder
maindir <- "path2benchmarkrnaseq/"
resp <- paste0(maindir, "roastgsa_on_benchmarking_data_rnaseq/")
dir.create(resp)


####### load data
## 1,2,6,7,8,9,10,12,14,15,16,17,20 - only GROUP
## 3,4,5,11,13,18,19, GROUP, BLOCK
tcga <- loadEData("tcga", nr.datasets=24)#,cache=FALSE)
names(tcga)

####### vector to save results
ares <- list()
RESULTS  <- list()

Ns  <- list()

####### Rotations and cores
nrot <- 1000
mccores <- 10

####### MalaCards scores
mala.kegg <- readRDS("KEGG.rds")
d2d.map <- readDataId2diseaseCodeMap("GseId2Disease.txt")
data(kegg.hs)

####### Execute roastgsa for every dataset
for(whdata in 1:length(tcga)){
    print(whdata)

    ysel <- assays(tcga[[whdata]])$expr
    name1 <- names(tcga)[whdata]
    name1%in%names(mala.kegg)


    dds1 <- DESeqDataSetFromMatrix(countData=ysel,colData=colData(tcga[[whdata]]),design= ~ GROUP)
    dds1 <- estimateSizeFactors(dds1)

    ynorm <- assays(vst(dds1))[[1]]
    ynorm <- ynorm[apply(ysel,1,mean)>10,]
    N  <- ncol(ysel)
    Ns[whdata]  <- N


    cnames <- c("BLOCK","GROUP")
    covar <- data.frame(colData(tcga[[whdata]])[,cnames,drop=F])
    covar$BLOCK <- as.factor(covar$BLOCK)
    covar$GROUP <- as.factor(covar$GROUP)
    colnames(covar) <- cnames
    index <- lapply(kegg.gs,function(x) rownames(ynorm)[which(rownames(ynorm)%in%x)])
    index  <- index[sapply(index,length)>5]

    design <- model.matrix( as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), data = covar)
    terms <- colnames(design)
    contrast <- which(colnames(design) == terms[length(terms)])


    ares[["absmean"]] <- roastgsa(ynorm, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                              index = index, nrot = nrot, mccores = mccores, set.statistic = "absmean",
                              self.contained = FALSE)


    ares[["maxmean"]] <- roastgsa(ynorm, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                  index = index, nrot = nrot, mccores = mccores, set.statistic = "maxmean",
                              self.contained = FALSE)


    ares[["mean"]] <- roastgsa(ynorm, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                  index = index, nrot = nrot, mccores = mccores,  set.statistic = "mean",
                              self.contained = FALSE)


    ares[["med"]] <- roastgsa(ynorm, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                              index = index, nrot = nrot, mccores = mccores, set.statistic = "median",
                              self.contained = FALSE)

    ares[["ksmean"]] <- roastgsa(ynorm, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                 index = index, nrot = nrot, mccores = mccores, set.statistic = "ksmean",
                             self.contained = FALSE)

    ares[["ksmax"]] <- roastgsa(ynorm, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                 index = index, nrot = nrot, mccores = mccores, set.statistic = "ksmax",
                             self.contained = FALSE)

    ares[["meanrank"]] <- roastgsa(ynorm, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                  index = index, nrot = nrot, mccores = mccores, set.statistic = "mean.rank",
                              self.contained = FALSE )

    ## measures M1 and M2
    RESULTS[[whdata]]  <- lapply(ares, function(x){
    a1  <- x$res
    a1$rank  <- rank(abs(a1[,which(regexpr("pval",names(x$res))>0)[1]]))
    rownames(a1)  <- substring(rownames(a1), 1, 8)
    M1  <- sum((1-a1[rownames(mala.kegg[[name1]]),]$rank/nrow(a1)) *
        mala.kegg[[name1]]$REL.SCORE, na.rm=TRUE)
    a2 <- a1
    pv <- sapply(1:1000,function(ds){
        rownames(a2) <- sample(rownames(a1))
        sum((1-a2[rownames(mala.kegg[[name1]]),]$rank/nrow(a1)) *
                mala.kegg[[name1]]$REL.SCORE, na.rm=TRUE)
    })
    M1.pv <- mean(M1<pv)
    whe  <- which(a1[rownames(mala.kegg[[name1]]),]$rank < 50)
    M2  <-  sum((1-a1[rownames(mala.kegg[[name1]]),]$rank/nrow(a1))[whe] *
        mala.kegg[[name1]]$REL.SCORE[whe], na.rm=TRUE)

    a2 <- a1
    pv <- sapply(1:1000,function(ds){
        rownames(a2) <- sample(rownames(a1))
        whe  <- which(a2[rownames(mala.kegg[[name1]]),]$rank < 50)
        sum((1-a2[rownames(mala.kegg[[name1]]),]$rank/nrow(a2))[whe] *
                mala.kegg[[name1]]$REL.SCORE[whe], na.rm=TRUE)
    })
    M2.pv <- mean(M2<pv)

    c(M1,M2, M1.pv, M2.pv)
   })
}

saveRDS(RESULTS, file = paste0(resp,"/RESULTS_rnaseq.RDS"))


#######  barplots
RESULTS <- readRDS(paste0(resp,"/RESULTS_rnaseq.RDS"))


pvs <- sapply(seq(0.05,0.5,by=0.05), function(K) apply(sapply(RESULTS,function(x) sapply(x,"[",3))<K,1,sum))
dff1 <- data.frame(sign.cases = as.numeric(pvs)+runif(length(as.numeric(pvs)),-0.2,0.2), type = rep(rownames(pvs),ncol(pvs)),   threshold = rep(seq(0.05,0.5,by=0.05),each=nrow(pvs)))

cols <- c("blue","green","orange","darkgreen","red","purple","pink")
g <- ggplot(dff1, aes(x=threshold,y=sign.cases, colour= type)) + geom_line(lwd=2) + geom_point(size=3)+theme_classic() + scale_colour_manual(values=cols) +
    geom_vline(xintercept=c(0.05,0.1,0.15,0.25),lty=2, colour="grey")+scale_y_continuous(name ="pvalue < threshold", breaks =c(1:11))
ggsave(g, file = paste0(resp, "sog_rnaseq.pdf"), width =10,height=7)


sapply(RESULTS,function(x) sapply(x,"[",3))

sum1  <- sapply(1:7, function(a) apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",1)),2,rank),1,function(x) sum(x==a)))
colnames(sum1)  <- paste0("rank_",1:7)
M1  <- apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",1)),2,rank),1,mean)

pdf(paste0(resp,"M1_barplot_rnaseq.pdf"), width =8, height =8)
layout(cbind(1:8,9:16),widths = c(1,0.1), height = c(rep(1,7),0.5))
par(mar=c(0,1.5,1,0))
for(k in 1:7) {
   bp  <-  barplot(t(sum1)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,14), col =k, ylab ="")
    title(ylab=rownames(sum1)[k], line=0, cex.lab=1.7 )
}
par(mar=c(0,1.5,1,0))
bp  <-  barplot(t(sum1)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,14), ylab ="",col="white", border=NA)
text(bp[,1], rep(6,7),paste0("rank_",1:7), cex = 1.5)

par(mar=c(0,0,0,0))
for(k in 1:7){
    plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
    text(0,-0.5,round(M1,2)[k], cex=1.7, col =1)
}
plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
text(0,-0.3,"avg", cex=1.7)
dev.off()

sum2  <- sapply(1:7, function(a) apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",2)),2,rank),1,function(x) sum(x==a)))
colnames(sum2)  <- paste0("rank_",1:7)
M2 <- apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",2)),2,rank),1,mean)

pdf(paste0(resp,"M2_barplot_rnaseq.pdf"), width =8, height =8)
layout(cbind(1:8,9:16),widths = c(1,0.1), height = c(rep(1,7),0.5))
par(mar=c(0,1.5,1,0))
for(k in 1:7) {
   bp  <-  barplot(t(sum2)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,13), col =k, ylab ="")
    title(ylab=rownames(sum2)[k], line=0, cex.lab=1.7 )
}
par(mar=c(0,1.5,1,0))
bp  <-  barplot(t(sum2)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,13), ylab ="",col="white", border=NA)
text(bp[,1], rep(6,7),paste0("rank_",1:7), cex = 1.5)

par(mar=c(0,0,0,0))
for(k in 1:7){
    plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
    text(0,-0.5,round(M2,2)[k], cex=1.7, col =1)
}
plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
text(0,-0.3,"avg", cex=1.7)
dev.off()


###### heatmaps
library(ComplexHeatmap)
library(circlize)

dad <-(apply(sapply(RESULTS,function(x) -sapply(x,"[",2)),2,rank)) -4
colnames(dad) <- rep("",ncol(dad))

pdf(paste0(resp,"heatmapwithalldata_M2_rnaseq.pdf"), height =5, width=10)
Heatmap((dad),  col = colorRamp2(c(3, 0, -3), c("blue", "white", "red")),
        heatmap_legend_param = list(title = "centred rank") )
dev.off()


dad <-(apply(sapply(RESULTS,function(x) -sapply(x,"[",1)),2,rank)) -4
colnames(dad) <- rep("",ncol(dad))

pdf(paste0(resp,"heatmapwithalldata_M1_rnaseq.pdf"), height =5, width=10)
Heatmap((dad),  col = colorRamp2(c(3, 0, -3), c("blue", "white", "red")),
                heatmap_legend_param = list(title = "centred rank") )
dev.off()



