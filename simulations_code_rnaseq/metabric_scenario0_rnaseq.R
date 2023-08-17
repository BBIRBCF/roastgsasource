#################################################################################################
##############  Scenario 0. Simulate new data from GTEX - breast: no effect  #####################
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

index2 <- readRDS(file = paste0(resp, "index2_sc0.RDS"))

####### Simulations
k <- 1
power1  <- list()
power1[[1]]  <- list()
power1[[2]]  <- list()
power1[[3]]  <- list()

for(k in c(1:length(Nselall))){

    ale <- mclapply(1:Nite, function(ite){
        set.seed(213 + ite)
        t1 <- thin_2group(as.matrix(counts.breast),group_prop = 0.5, signal_params = list(mean = 0.2, sd = 1),prop_null = 1)
        samps <- c(sample(which(t1$designmat==1),Nselall[k]/2),sample(which(t1$designmat==0),Nselall[k]/2))
        c1 <- cor(t(ynormall[index2[[1]],]))
        (avg.cor <- mean(c1[lower.tri(c1)]))

        ysim <- t1$mat[,samps]
        N  <- ncol(ysim)

        group <- t1$designmat[samps]

        group <- data.frame(group=as.factor(group))
        rownames(group) <- paste0("s",1:nrow(group))
        rownames(ysim) <- rownames(counts.breast)
        colnames(ysim) <- rownames(group)

        dds1 <- DESeqDataSetFromMatrix(countData=ysim,colData=group,design= ~  group)
        dds1 <- DESeq(dds1)
        if(k==5)
            ynorm <- assays(vst(dds1))[[1]]
        else
            ynorm <- assays(rlog(dds1))[[1]]

        (CO.median <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                              index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                              set.statistic = "median")$res[names(index2),]$pval)


        SC.median  <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                               index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                               set.statistic = "median")$res[names(index2),]$pval


        SC.absmean  <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                                index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                                set.statistic = "absmean")$res[names(index2),]$pval

        CO.absmean  <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                                index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                                set.statistic = "absmean")$res[names(index2),]$pval

        CO.mean.rank  <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                                  index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                                  set.statistic = "mean.rank")$res[names(index2),]$pval.diff

        (SC.maxmean  <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                                index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                                set.statistic = "maxmean")$res[names(index2),]$pval)


        CO.maxmean <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                               index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                               set.statistic = "maxmean")$res[names(index2),]$pval

        CO.ksmean <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                            index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                            set.statistic = "ksmean")$res[names(index2),]$pval

        CO.ksmax  <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                             index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                             set.statistic = "ksmax")$res[names(index2),]$pval

        SC.mean <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                            index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                            set.statistic = "mean")$res[names(index2),]$pval
        CO.mean  <- roastgsa((ynorm), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(group[,1])),
                             index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                             set.statistic = "mean")$res[names(index2),]$pval

        list(SC.mean, SC.maxmean, SC.median, SC.absmean,
             CO.mean, CO.maxmean, CO.median, CO.absmean, CO.mean.rank, CO.ksmean, CO.ksmax)

    }, mc.cores = mc.cores)


    power1[[1]][[k]] <- lapply(1:length(ale[[1]]),function(o) sapply(ale, function(x) x[[o]][1]))
    power1[[2]][[k]] <- lapply(1:length(ale[[1]]),function(o) sapply(ale, function(x) x[[o]][2]))
    power1[[3]][[k]] <- lapply(1:length(ale[[1]]),function(o) sapply(ale, function(x) x[[o]][3]))
   save(power1, avgcor2t, file= paste0(resp, "scenario0_empiricalsize_rnaseq.Rdata"))
}

###### Write recovery rates into table
load(paste0(resp, "scenario0_empiricalsize_rnaseq.Rdata"))
nm  <- c("SC.mean", "SC.maxmean",  "SC.median", "SC.absmean",
         "CO.mean", "CO.maxmean", "CO.median", "CO.absmean", "CO.mean.rank", "CO.ksmean", "CO.ksmax")

cat("\n", file = paste0(resp, "scenario0_table.txt"), append = FALSE)

for(K in c(1,2,3)){
    print(K)
    power2 <- power1[[K]]
    avg.tab <-  do.call(cbind,lapply(1:5,function(o) sapply(1:11, function(jk)
        mean(power2[[o]][[jk]]<0.05))))
    rownames(avg.tab) <- nm
    colnames(avg.tab) <- paste0("N.",Nselall)
    tab1  <- avg.tab[nm,] *1000
    rownames(tab1)  <- nm

    for(o in 1:nrow(tab1))
    {
        if(o==5) cat(paste0("A",K," & ",rownames(tab1)[o], "& ", paste0(tab1[o,], collapse =" & "), '\\','\\'), "\n", file = paste0(resp, "scenario0_table.txt"), append = TRUE)
        else cat(paste0("& ",rownames(tab1)[o], "& ", paste0(tab1[o,], collapse =" & "), '\\','\\'), "\n", file = paste0(resp, "scenario0_table.txt"), append = TRUE)
    }
    cat("\\hline ", "\n", file = paste0(resp, "scenario0_table.txt"), append = TRUE)
}

