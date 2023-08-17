#################################################################################################
##############  Scenario 1. Simulate new data from GTEX - breast: same effect- low cor   ########
#################################################################################################
## R version 4.1.3 (2023-05-10)
library(roastgsa)
library(seqgendiff)
library(limma)
library(DESeq2)

notdonebefor  <- FALSE
Nite  <- 1000
nrot  <- 500
mc.cores  <- 5

maindir <- "path2maindir/"
dir.create(resp <- paste0(maindir,"/scenarios_recoveryrates_rnaseq/"))

counts.breast <- readRDS(file = paste0(res, "counts.breast.RDS"))
pd.breast <- readRDS(file = paste0(res, "pd.breast.RDS"))

####### Number of samples and testing sets
Nselall <- c(6, 10, 20, 30, 100)

cda <- as.data.frame(rep(1, ncol(counts.breast)))
rownames(cda) <- colnames(counts.breast)

dds1 <- DESeqDataSetFromMatrix(countData=counts.breast,colData=cda,design= ~  1)
dds1 <- estimateSizeFactors(dds1)
ynormall <- assays(vst(dds1))[[1]]
rownames(ynormall) <- rownames(counts.breast)
basemean <- apply(counts.breast,1,mean)
ynormall2 <- ynormall[basemean>10,]

####### Simulations
k <- 1
avgcor2t <- numeric()
power1  <- list()
power1[[1]]  <- list()
power1[[2]]  <- list()
power1[[3]]  <- list()

index <- list(i1= sample(rownames(ynormall2),40), i2=sample(rownames(ynormall2),40))
for(k in 1:length(Nselall)){

    ale <- mclapply(1:Nite, function(ite){
        print(ite)
        set.seed(213 + ite)

        s1 <- sample(rownames(ynormall2),30)

        t1 <- thin_2group(as.matrix(counts.breast),group_prop = 0.5, signal_params = list(mean = 0.15, sd = 0.001),prop_null = 0,alpha=1)
        samps <- c(sample(which(t1$designmat==1),Nselall[k]/2),sample(which(t1$designmat==0),Nselall[k]/2))

        bm1 <- apply(counts.breast[,samps],1,mean)

        index2 <- index
        index2[[1]] <- s1

        c1 <- cor(t(ynormall[index2[[1]],]))
        avg.cor <- mean(c1[lower.tri(c1)])

        aux <- t1$mat
        rownames(aux) <- rownames(counts.breast)
        ysim <- counts.breast[,samps]
        ysim[index2[[1]],] <- aux[index2[[1]],samps]
        N  <- ncol(ysim)
        basemean <- apply(ysim,1,mean)[index2[[1]]]
        group <- t1$designmat[samps]

        group <- data.frame(group=as.factor(group))
        rownames(group) <- paste0("s",1:nrow(group))
        rownames(ysim) <- rownames(counts.breast)
        colnames(ysim) <- rownames(group)

        dds1 <- DESeqDataSetFromMatrix(countData=ysim,colData=group,design= ~  group)
        dds1 <- estimateSizeFactors(dds1)
        ynorm <- assays(vst(dds1))[[1]]
        ynorm2 <- ynorm[apply(ysim,1,mean)>10,]
        index2 <- lapply(index2, function(x) x[x%in%rownames(ynorm2)])

        (CO.median <- roastgsa((ynorm2), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                              index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                              set.statistic = "median")$res[names(index2),]$pval)

        CO.absmean  <- roastgsa((ynorm2), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                                index =index2 , nrot = nrot, mccores = 1, self.contained = FALSE,
                                set.statistic = "absmean")$res[names(index2),]$pval

        CO.mean.rank  <- roastgsa((ynorm2), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                                  index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                                  set.statistic = "mean.rank")$res[names(index2),]$pval.diff

        CO.ksmean <- roastgsa((ynorm2), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                            index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                            set.statistic = "ksmean")$res[names(index2),]$pval

        CO.ksmax  <- roastgsa((ynorm2), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                             index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                             set.statistic = "ksmax")$res[names(index2),]$pval

        CO.mean  <- roastgsa((ynorm2), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                             index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                             set.statistic = "mean")$res[names(index2),]$pval

        CO.maxmean <- roastgsa((ynorm2), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                               index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                               set.statistic = "maxmean")$res[names(index2),]$pval

        list(CO.mean, CO.maxmean, CO.median, CO.absmean, CO.mean.rank, CO.ksmean, CO.ksmax)

    }, mc.cores = mc.cores)

    power1[[1]][[k]] <- lapply(1:length(ale[[1]]),function(o) sapply(ale, function(x) x[[o]][1]))
    save(power1, basemean, avgcor2t, file= paste0(resp, "scenario1_lowavgcor_rnaseq.Rdata"))
}




load(paste0(resp, "scenario1_lowavgcor_rnaseq.Rdata"))

nm  <- c("CO.mean", "CO.maxmean", "CO.median", "CO.absmean", "CO.mean.rank", "CO.ksmean", "CO.ksmax")

cat("\n", file = paste0(resp, "scenario1_lowavgcor.txt"), append = FALSE)

for(K in c(1)){
    print(K)
    power2 <- power1[[K]]
    avg.tab <-  do.call(cbind,lapply(1:5,function(o) sapply(1:7, function(jk)
        mean(power2[[o]][[jk]]<0.05))))
    rownames(avg.tab) <- nm
    colnames(avg.tab) <- paste0("N.",Nselall[1:5])
    tab1  <- avg.tab[nm,] *1000
    rownames(tab1)  <- nm

    for(o in 1:nrow(tab1))
    {
        if(o==5) cat(paste0("A",K," & ",rownames(tab1)[o], "& ", paste0(tab1[o,], collapse =" & "), '\\','\\'), "\n", file = paste0(resp, "scenario1_lowavgcor.txt"), append = TRUE)
        else cat(paste0("& ",rownames(tab1)[o], "& ", paste0(tab1[o,], collapse =" & "), '\\','\\'), "\n", file = paste0(resp, "scenario1_lowavgcor.txt"), append = TRUE)
    }
    cat("\\hline ", "\n", file = paste0(resp, "scenario1_lowavgcor.txt"), append = TRUE)
}
