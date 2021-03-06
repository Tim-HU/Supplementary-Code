###################
## Integrating snDrop-seq and scTHS-seq data
## @author: Jean Fan
## @email: jeanfan@fas.harvard.edu
###################

require(dendextend)
library("Rcpp", lib.loc="/usr/local/lib/R/site-library")
library(WGCNA)
library("dbscan", lib.loc="/usr/local/lib/R/site-library")
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(Cairo)
library(parallel)
library("pagoda2", lib.loc="/home/barkasn/R/x86_64-pc-linux-gnu-library/3.3/") 
source("/home/pkharchenko/m/p2/fastclust.r")
source("/home/pkharchenko/m/p2/sparseclust.r")
source("/home/pkharchenko/m/p2/schromoclust.r")

########################## Occ integrative analysis

## snDrop-seq results
load("fin.blue.06152017.RData")
DS.r <- Occ
DS.annot <- factor(fin.blue[rownames(DS.r$counts)])
table(DS.annot)
lvec <- colSumByFac(DS.r$misc[['rawCounts']][names(DS.annot),],DS.annot)[-1,] + 1
lvec <- t(lvec/pmax(1,rowSums(lvec)))
colnames(lvec) <- levels(DS.annot)
rownames(lvec) <- colnames(DS.r$misc[['rawCounts']])
ld <- jsDist(lvec); colnames(ld) <- rownames(ld) <- colnames(lvec)
hctree <- stats::hclust(as.dist(ld),method='ward.D')
plot(hctree)
DS.tree <- hctree

## scTHS-seq results
load('~/Projects/Kun_Epigenetics/R5/hR10_vc_comb_fixed_withopcend.RData')
THS.r <- r
THS.annot <- annot
lvec <- colSumByFac(THS.r$misc[['rawCounts']],THS.annot)[-1,] + 1
lvec <- t(lvec/pmax(1,rowSums(lvec)))
colnames(lvec) <- levels(THS.annot)
rownames(lvec) <- colnames(THS.r$misc[['rawCounts']])
ld <- jsDist(lvec); colnames(ld) <- rownames(ld) <- colnames(lvec)
hctree <- stats::hclust(as.dist(ld),method='ward.D')
plot(hctree)
THS.tree <- hctree

par(mfrow=c(2,2), mar=rep(5,4))
DS.r$plotEmbedding(type='PCA',show.legend=F,groups=DS.annot,mark.clusters=T,mark.cluster.cex=1, main='snDrop-seq',embeddingType='tSNE', alpha=0.1)
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=THS.annot,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1)
plot(DS.tree, axes=FALSE)
plot(THS.tree, axes=FALSE)

############################################### predict
## Annotate peaks with genes
peaks <- colnames(THS.r$counts)
peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
peak.df$mid <- (as.numeric(peak.df$end)+as.numeric(peak.df$start))/2
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb)
library(ChIPseeker)
peak.gr <- with(peak.df, GRanges(as.character(chr), IRanges(as.numeric(as.character(start)), as.numeric(as.character(end)))))
p2g <- annotatePeak(peak.gr, TxDb=txdb)

p2g.detail <- cbind(geneId=p2g@anno$geneId, distanceToTSS=p2g@anno$distanceToTSS, p2g@detailGenomicAnnotation)
rownames(p2g.detail) <- peaks

require(annotate)
library(org.Hs.eg.db)
p2g.detail$symbol <- getSYMBOL(p2g.detail$geneId, data='org.Hs.eg')
head(p2g.detail)

############################################# Training model
## use Ast vs. Oli to learn features
g1 <- 'Ast'
g2 <- 'Oli'
test <- c(g1, g2)
direct.annot <- as.character(DS.annot)
direct.annot[!(direct.annot %in% test)] <- NA
table(direct.annot)
direct.annot <- factor(direct.annot)
names(direct.annot) <- names(DS.annot)
DS.r$plotEmbedding(type='PCA',show.legend=F,groups=direct.annot,mark.clusters=T,mark.cluster.cex=1, main='snDrop-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)
direct.dg <- DS.r$getDifferentialGenes(upregulated.only = TRUE, groups=direct.annot, z.threshold=0.01)

topAstGenes <- direct.dg[[g1]]
topAstGenes <- rownames(topAstGenes)[topAstGenes$Z>1.96]
topAstGenes <- intersect(topAstGenes, p2g.detail$symbol)
length(topAstGenes)
topOliGenes <- direct.dg[[g2]]
topOliGenes <- rownames(topOliGenes)[topOliGenes$Z>1.96]
topOliGenes <- intersect(topOliGenes, p2g.detail$symbol)
length(topOliGenes)

DS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=Matrix::rowSums(DS.r$counts[,topAstGenes]),alpha=0.1,main="TopAstPeaks",embeddingType='tSNE')
DS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=Matrix::rowSums(DS.r$counts[,topOliGenes]),alpha=0.1,main="TopOliPeaks",embeddingType='tSNE')

## select sites that are near diff genes
topAstPeaks <- rownames(p2g.detail)[p2g.detail$symbol %in% topAstGenes]
topOliPeaks <- rownames(p2g.detail)[p2g.detail$symbol %in% topOliGenes]
THS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=Matrix::rowSums(THS.r$counts[,topAstPeaks]),alpha=0.1,main="TopAstPeaks",embeddingType='tSNE')
THS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=Matrix::rowSums(THS.r$counts[,topOliPeaks]),alpha=0.1,main="TopOliPeaks",embeddingType='tSNE')

## get real diff sites
direct.annot <- as.character(THS.annot)
direct.annot[!(direct.annot %in% test)] <- NA
table(direct.annot)
direct.annot <- factor(direct.annot)
names(direct.annot) <- names(THS.annot)
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=direct.annot,mark.clusters=T,mark.cluster.cex=1, main='snDrop-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)
direct.da <- THS.r$getDifferentialGenes(upregulated.only = TRUE, groups=direct.annot, z.threshold=0.01)

## train classifier to predict diff peaks
sdat <- p2g.detail

table(rownames(direct.dg[[1]]) %in% rownames(direct.dg[[2]]))  ## mutually exclusive
sdat.dg <- rbind(direct.dg[[g1]], direct.dg[[g2]])
#colnames(sdat.dg) <- c(paste0('dgg1.', colnames(direct.dg[[g1]])), paste0('dgg2.', colnames(direct.dg[[g2]])))
head(sdat.dg)

table(names(direct.da[[1]]) %in% names(direct.da[[2]]))  ## mutually exclusive
vi <- names(direct.da[[1]]) %in% names(direct.da[[2]])
direct.da[[1]] <- direct.da[[1]][!vi]
direct.da[[2]] <- direct.da[[2]][!vi]
sdat.da <- c(direct.da[[g1]], direct.da[[g2]])
head(sdat.da)

sdat.all <- cbind(sdat.da, p2g.detail[names(sdat.da),])
table(sdat.all$symbol %in% rownames(sdat.dg))
sdat.all <- cbind(sdat.all, sdat.dg[sdat.all$symbol,])

head(sdat.all)
sdat.all$diff <- ifelse(sdat.all$sdat.da>1.96, 'diff','nodiff')
table(sdat.all$diff)
head(sdat.all)

vi <- !is.na(sdat.all$Z)
table(vi)

sdat.filter <- sdat.all[vi,]
table(sdat.filter$diff)
sdat.filter$sdat.da <- NULL
sdat.filter$geneId <- NULL
sdat.filter$symbol <- NULL
head(sdat.filter, n=10)

## Train
dim(sdat.filter)
sdat.sub <- sdat.filter[sample(1:nrow(sdat.filter),6e3),]
table(sdat.sub$diff)
head(sdat.sub)

## for a single site, predict if it will be differential based on gene info
library(caret)
fitControl <- trainControl( method = "cv", number = 10, classProbs=TRUE, savePredictions=T,summaryFunction=twoClassSummary)
epi.gbmFit <- train(diff ~ ., data = sdat.sub,  method = "gbm",  trControl = fitControl, metric='ROC', verbose = FALSE)
#epi.gbmFit <- train(diff ~ ., data = sdat.sub,  method = "svmRadial",  trControl = fitControl, metric='ROC', verbose = FALSE)
names(epi.gbmFit)
head(epi.gbmFit$pred$obs)
head(epi.gbmFit$pred$diff)

save(epi.gbmFit, p2g.detail, file='epi.gbmFit.ast.oli.RData')

library(pROC)
par(mfrow=c(2,2))
plot.roc(epi.gbmFit$pred$obs, epi.gbmFit$pred$diff,percent=T,ci=T,print.auc=T)

## Look at feature weights
importance <- varImp(epi.gbmFit, scale=FALSE)
print(importance)

## Jointly
pdf <- sdat.all[sdat.all$symbol %in% topAstGenes,]
pdf$sdat.da <- NULL
pdf$geneId <- NULL
pdf$symbol <- NULL
head(pdf)
pre <- predict(epi.gbmFit,pdf,type='prob')
#x <- pre$diff - pre$nodiff; names(x) <- rownames(pdf)
x <- pre$diff; names(x) <- rownames(pdf)
range(x)
pAst <- x

pdf <- sdat.all[sdat.all$symbol %in% topOliGenes,]
pdf$sdat.da <- NULL
pdf$geneId <- NULL
pdf$symbol <- NULL
head(pdf)
pre <- predict(epi.gbmFit,pdf,type='prob')
#x <- pre$diff - pre$nodiff; names(x) <- rownames(pdf)
x <- pre$diff; names(x) <- rownames(pdf)
range(x)
pOli <- x

range(pAst)
range(pOli)

t <- 0
dAst <- Matrix::colMeans(t(THS.r$counts[,names(pAst)[pAst>t]])*pAst[pAst>t])
dOli <- Matrix::colMeans(t(THS.r$counts[,names(pOli)[pOli>t]])*pOli[pOli>t])

pdscale <- rbind(dOli, dAst)
pdscale <- scale(pdscale, scale=FALSE)
#pdscale <- t(scale(t(pdscale), scale=FALSE))
pdscale[, names(THS.annot)[which(is.na(THS.annot))]] <- NA
pdscale[, names(direct.annot)[!(direct.annot %in% test)]] <- NA

## ROC curve
library(pROC)
foo <- pdscale[1,]-pdscale[2,]
f <- direct.annot[names(foo)]
rocobj <- plot.roc(f,foo,percent=T,ci=T,print.auc=T)

#par(mfrow=c(1,2))
THS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=pdscale[1,],alpha=0.1,main="dOli",embeddingType='tSNE')
THS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=pdscale[2,],alpha=0.1,main="dAst",embeddingType='tSNE')






#################### make into function
getSdat <- function(direct.dg, direct.annot) {
    print(table(rownames(direct.dg[[1]]) %in% rownames(direct.dg[[2]])))  ## mutually exclusive
    sdat.dg <- rbind(direct.dg[[1]], direct.dg[[2]])
    head(sdat.dg)

    sdat.all <- p2g.detail    
    table(sdat.all$symbol %in% rownames(sdat.dg))
    sdat.all <- cbind(sdat.all, sdat.dg[sdat.all$symbol,])

    pdf <- sdat.all[sdat.all$symbol %in% rownames(direct.dg[[1]]),]
    pdf$sdat.da <- NULL
    pdf$geneId <- NULL
    pdf$symbol <- NULL
    head(pdf)
    pre <- predict(epi.gbmFit,pdf,type='prob')
    x <- pre$diff; names(x) <- rownames(pdf)
    range(x)
    pAst <- x

    pdf <- sdat.all[sdat.all$symbol %in% rownames(direct.dg[[2]]),]
    pdf$sdat.da <- NULL
    pdf$geneId <- NULL
    pdf$symbol <- NULL
    head(pdf)
    pre <- predict(epi.gbmFit,pdf,type='prob')
    x <- pre$diff; names(x) <- rownames(pdf)
    range(x)
    pOli <- x
    
    t <- 0
    dAst <- Matrix::colMeans(t(THS.r$counts[,names(pAst)[pAst>t]])*pAst[pAst>t])
    dOli <- Matrix::colMeans(t(THS.r$counts[,names(pOli)[pOli>t]])*pOli[pOli>t])

    pdscale <- rbind(dAst, dOli)
    pdscale <- scale(pdscale, scale=FALSE)
    pdscale[, names(THS.annot)[which(is.na(THS.annot))]] <- NA
    
    return(pdscale)
}

getSdatInd <- function(direct.dg) {
    sdat.dg <- direct.dg
    head(sdat.dg)

    sdat.all <- p2g.detail    
    table(sdat.all$symbol %in% rownames(sdat.dg))
    sdat.all <- cbind(sdat.all, sdat.dg[sdat.all$symbol,])

    pdf <- sdat.all[sdat.all$symbol %in% rownames(direct.dg),]
    pdf$sdat.da <- NULL
    pdf$geneId <- NULL
    pdf$symbol <- NULL
    head(pdf)
    pre <- predict(epi.gbmFit,pdf,type='prob')
    x <- pre$diff; names(x) <- rownames(pdf)
    
    return(x)
}








##################### recursive classification
classify <- function(tree.init, cells.init, plot=FALSE, refine=1, p=0.9, t=0.9, num.cells=40, zs=1.28) {
    tryCatch({
    t1 <- tree.init[[1]]
    t2 <- tree.init[[2]]
    c1 <- labels(t1)
    c2 <- labels(t2)
    direct.cut <- c(rep(1, length(c1)), rep(2, length(c2))); names(direct.cut) <- c(c1,c2)
    
    ## new level group names
    newnames <- direct.cut
    namesmap <- sapply(unique(newnames), function(x) {
        paste0(names(newnames[newnames==x]), collapse=".")
    })
    names(namesmap) <- unique(newnames)
    namesfinal <- namesmap[newnames]
    names(namesfinal) <- names(newnames)
    table(namesfinal)

    ## map from DS to THS
    direct.annot <- factor(DS.annot, levels=levels(DS.annot), labels=namesfinal[levels(DS.annot)])
    direct.annot <- factor(direct.annot)
    table(direct.annot)
    direct.annot.all <- direct.annot
    ## sample equal numbers from subgroups
    set.seed(0)
    testcells <- lapply(names(namesfinal), function(on) {
        sample(names(na.omit(DS.annot[DS.annot==on])), num.cells, replace=FALSE)
    })
    names(testcells) <- names(namesfinal)
    testcells <- unlist(testcells)
    direct.annot[!(names(direct.annot) %in% testcells)] <- NA
    table(direct.annot)
    ## differentially expressed genes
    direct.dg <- DS.r$getDifferentialGenes(upregulated.only = TRUE, groups=direct.annot, z.threshold=zs)
    names(direct.dg)

    if(plot) {
        par(mfrow=c(1,1), mar=rep(5,4))
        c1 <- rowMeans(DS.r$counts[names(na.omit(DS.annot)),rownames(direct.dg[[1]])])
        c2 <- rowMeans(DS.r$counts[names(na.omit(DS.annot)),rownames(direct.dg[[2]])])
        col <- scale(c1) - scale(c2)
        col <- col[,1]
        col[col < -2] <- -2
        col[col > 2] <- 2
        col[!(names(col) %in% names(na.omit(direct.annot.all)))] <- 0
        DS.r$plotEmbedding(type='PCA',embedding='tSNE',show.legend=F,colors=col,alpha=0.1)
    }

    pd <- getSdat(direct.dg, direct.annot)
    rownames(pd) <- names(direct.dg)
    dim(pd)
    pd[is.nan(pd)] <- 0
    
    pdscale <- scale(pd)
    #pdscale <- pdscale[, cells.init]
    
    if(plot) {
        par(mfrow=c(1,1), mar=rep(5,4))
        c1 <- pdscale[1,][rownames(THS.r$counts)]
        c2 <- pdscale[2,][rownames(THS.r$counts)]
        col <- scale(c1)-scale(c2)
        col <- col[,1]
        #col <- c1 - c2
        col[!(names(col) %in% cells.init)] <- NA
        col[col < -2] <- -2
        col[col > 2] <- 2
        names(col) <- rownames(THS.r$counts)
        col[is.na(col)] <- 0
        THS.r$plotEmbedding(type='PCA',embedding='tSNE',show.legend=F,colors=col,alpha=0.1, main=rownames(pdscale)[1])
    }
    
    ## max vote
    vote <- factor(unlist(apply(pdscale, 2, function(x) {
        g <- names(x)[which(x==max(x, na.rm=TRUE))]
        if(max(x, na.rm=TRUE) < 0) { g <- NA }
        set.seed(0)
        sample(g, 1) ## randomly break ties
    })))
    vote <- na.omit(vote)
    vote <- vote[cells.init]
    print(table(vote))
    
    if(plot) {
        THS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,groups=vote,alpha=0.1,main="vote",embeddingType='tSNE')
    }
    

    for(i in 1:refine) {
        if(sum(table(vote)>0) != 2) {
            print('REFINEMENT FAILED')
            ## refinement caused issued or maybe other group not distinct
            #return(rep(NA, tree.init %>% nnodes))
            break;
            #return(list(g1=NA, g2=NA, t1=NA, t2=NA, rocs=NA, aucs=NA, vote=NA, names=paste0(unique(namesfinal), collapse=' vs ')))
        }

        dg <- THS.r$getDifferentialGenes(upregulated.only = TRUE, groups=vote, z.threshold=zs)
        lvec <- do.call(rbind, lapply(dg, function(g) {
            lv <- Matrix::rowSums(THS.r$counts[,names(g)])/length(g)/Matrix::rowSums(THS.r$counts)
            lv[is.na(vote)] <- NA
            lv
        }))
        lvec.scale <- t(scale(t(lvec)))

        if(plot) {
            par(mfrow=c(1,1), mar=rep(5,4))
            c1 <- lvec.scale[1,][rownames(THS.r$counts)]
            c2 <- lvec.scale[2,][rownames(THS.r$counts)]
            col <- scale(c1) - scale(c2)
            col <- col[,1]
            col[!(names(col) %in% cells.init)] <- NA
            col[col < -2] <- -2
            col[col > 2] <- 2
            THS.r$plotEmbedding(type='PCA',embedding='tSNE',show.legend=F,colors=col,alpha=0.1)
        }

        ## max vote
        pdscale <- lvec.scale
        vote <- factor(unlist(apply(pdscale, 2, function(x) {
            g <- names(x)[which(x==max(x, na.rm=TRUE))]
            if(max(x, na.rm=TRUE) < 0) { g <- NA }
            set.seed(0)
            sample(g, 1) ## randomly break ties
        })))
        vote <- na.omit(vote)
        vote <- vote[cells.init]
        if(plot) {
            THS.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,groups=vote,alpha=0.1,main="vote",embeddingType='tSNE')
        }
        
    }
    
    ## assess stability
    set.seed(1)
    ## sample evenly
    sample <- unlist(lapply(levels(vote), function(v) {
        gv <- names(vote[vote==v])
        sample(gv, length(gv)*p)
    }))
    #sample <- sample(names(vote), length(vote)*p)
    samplenot <- setdiff(names(vote), sample)
    
    if(sum(table(vote[sample])>0) <= 1) {
        print('STABILITY ASSESSMENT FAILED')
        ## refinement caused issued or maybe other group not distinct
        return(rep(NA, tree.init %>% nnodes))
        #return(list(g1=NA, g2=NA, t1=NA, t2=NA, rocs=NA, aucs=NA, vote=NA, names=paste0(unique(namesfinal), collapse=' vs ')))
    }

    dg <- THS.r$getDifferentialGenes(upregulated.only = TRUE, groups=vote[sample], z.threshold=zs)
    lvec <- do.call(rbind, lapply(dg, function(g) {
        lv <- Matrix::rowSums(THS.r$counts[,names(g)])/length(g)/Matrix::rowSums(THS.r$counts)
        lv[is.na(vote)] <- NA
        lv
    }))
    lvec <- lvec[, cells.init]    
    lvec.scale <- t(scale(t(lvec)))
    #lvec.scale[is.nan(lvec.scale)] <- NA
    #lvec.scale[is.infinite(lvec.scale)] <- NA
    
    if(plot) {
        par(mfrow=c(1,1), mar=rep(5,4))
        c1 <- lvec.scale[1,][rownames(THS.r$counts)]
        c2 <- lvec.scale[2,][rownames(THS.r$counts)]
        col <- scale(c1) - scale(c2)
        col <- col[,1]
        col[!(names(col) %in% cells.init)] <- NA
        col[col < -2] <- -2
        col[col > 2] <- 2
        THS.r$plotEmbedding(type='PCA',embedding='tSNE',show.legend=F,colors=col,alpha=0.1)
    }

    mat <- data.frame(t(lvec.scale))
    mat[, 'annot'] <- vote[rownames(mat)]
    head(mat)
    mat1 <- mat[sample,]
    mat1 <- na.omit(mat1)
    mat2 <- mat[samplenot,]
    mat2 <- na.omit(mat2)

    ## train on same cells used in diff exp
    require(caret)
    fitControl <- trainControl( method = "cv", number = 10, classProbs=TRUE, savePredictions=T)
    z1 <- train(annot ~ ., mat1, method = 'gbm', trControl = fitControl, verbose=FALSE)
    z2 <- predict(z1, newdata=mat2, type="prob")
    rownames(z2) <- rownames(mat2)
    ## assess performance based on cells not used for diff exp
    require(pROC)
    x <- levels(vote)[1]
    if(plot) {        
        rocs <- roc(mat2$annot==x, z2[,x], plot=TRUE, print.auc=TRUE, ci=TRUE)
    }
    else {
        rocs <- roc(mat2$annot==x, z2[,x], plot=FALSE)
    }
    aucs <- rocs$auc
    print(aucs)
    
    g1 <- na.omit(names(vote)[as.character(vote)==namesmap[1]])
    g2 <- na.omit(names(vote)[as.character(vote)==namesmap[2]])

    return(list(g1=g1, g2=g2, t1=t1, t2=t2, rocs=rocs, aucs=aucs, vote=vote, lvec=lvec, names=paste0(unique(namesfinal), collapse=' vs ')))
    }, error = function(e) {
        return(rep(NA, tree.init %>% nnodes))
    })
}

pdf('tree_example.pdf')
par(mfrow=c(1,1), mar=rep(5,4))
plot(as.dendrogram(DS.tree))
dev.off()

pdf('snDropseq2scTHSseq_example.pdf', width=3, height=3)

tree.init <- as.dendrogram(DS.tree)
cells.init <- names(na.omit(THS.annot))
results1 <- classify(tree.init, cells.init, plot=FALSE, refine=5)
test <- results1$vote
print(table(test))
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.2 <- classify(results1$t2, results1$g2, plot=FALSE, refine=5)
test <- results2.2$vote
print(table(test))
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

## results2.1 <- classify(results1$t1, results1$g1, plot=TRUE, refine=5)
## test <- results2.1$vote
## THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.2.1 <- classify(results2.2$t1, results2.2$g1, plot=TRUE, refine=5, zs=0.1)
#results2.2.1 <- classify(results2.2$t1, ex.cells, plot=TRUE, refine=5)
test <- results2.2.1$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.2.2 <- classify(results2.2$t2, results2.2$g2, plot=TRUE, refine=5, t=0.4)
#pdf('snDropseq2scTHSseq_example_in.pdf', width=3, height=3)
#results2.2.2 <- classify(results2.2$t2, in.cells, plot=TRUE, t=0.4, zs=3)
test <- results2.2.2$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.2.2.1 <- classify(results2.2.2$t1, results2.2.2$g1, plot=TRUE, t=0.4, zs=3)
test <- results2.2.2.1$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.2.2.2 <- classify(results2.2.2$t2, results2.2.2$g2, plot=TRUE, t=0.4)
test <- results2.2.2.2$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)
#dev.off()

dev.off()


################## Manual run

tree.init <- as.dendrogram(DS.tree)
cells.init <- names(na.omit(THS.annot))
results1 <- classify(tree.init, cells.init, plot=FALSE, refine=5, t=0.4, zs=1.28)
test <- results1$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.1 <- classify(results1$t1, results1$g1, plot=FALSE, refine=5, t=0.9, zs=1.28)
test <- results2.1$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.1.2 <- classify(results2.1$t2, results2.1$g2, plot=TRUE, refine=5, t=0.9, zs=1.96)
test <- results2.1.2$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.1.1 <- classify(results2.1$t1, results2.1$g1, plot=FALSE, refine=5, t=0.4, zs=3)
test <- results2.1.1$vote
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=test,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

results2.2.1.1 <- classify(results2.2.1$t1, results2.2.1$g1, plot=TRUE)
results2.2.1.2 <- classify(results2.2.1$t2, results2.2.1$g2, plot=TRUE)



############### Recursive algorithm
classify.recur <- function(tree.init, cells.init) {
    print(labels(tree.init))
    refine = 5
    zs = 3
    if(sum(grepl('Ast|Oli|End|Per|Mic|OPC', tree.init %>% labels)) > 2) {
       zs = 1.28
    } 
    if(length(labels(tree.init)) > 1) {
        results <- classify(tree.init, cells.init, refine=refine, zs=zs)
        if(sum(!is.na(results))==0) {
            refine = 5
            results <- classify(tree.init, cells.init, refine=refine, zs=1.28)
        }
        #results <- classify(tree.init, cells.init)
        all.results <<- c(all.results, list(results))

        ## recur
        if(!is.na(results)) {
            classify.recur(results$t1, results$g1)
            classify.recur(results$t2, results$g2)
        } 
        
    } else {
        ## is leaf, just store NA
        all.results <<- c(all.results, NA)
    }
}

all.results <- list()
tree.init <- as.dendrogram(DS.tree)
cells.init <- names(na.omit(THS.annot))
classify.recur(tree.init, cells.init)

## visualize as tree
bw <- unlist(lapply(all.results, function(x) {
    if(!is.na(x)) {
        y <- x$aucs
        names(y) <- x$name
        y
    } else {
        x
    }
}))
bw
bw[is.na(bw)] <- 0
length(bw)
tree.init %>% nnodes

## scale colors
bw[bw < 0.5] <- 0.5
bwn <- round(bw*9)+1
colcols <- colorRampPalette(c('grey', 'yellow', 'red'))(10)[bwn]

library(dendextend)
par(mfrow=c(1,1), mar=rep(5,4))
#pdf('final_tree_probs2.pdf', width=6, height=6)
tree.init %>% set("nodes_pch", 19) %>% set("nodes_col", colcols) %>% plot()
#dev.off()

#pdf('final_cluster_probs_all2.pdf')
bw <- do.call(cbind, lapply(all.results, function(x) {
    classification <- rep(NA, length(all.cells))
    names(classification) <- all.cells
    if(!is.na(x)) {
        y <- x$auc
        if(y > 0.7) {
            classification[names(x$vote)] <- as.character(x$vote)
            THS.r$plotEmbedding(type='PCA',show.legend=F,groups=classification,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1)            
        }
    }
    return(classification)
}))
#dev.off()
head(bw)

final <- apply(bw, 1, function(x) {
    if(sum(!is.na(x))>0) {
        y <- na.omit(x)
        ## shortest name is finest classification?    
        z <- sapply(y, function(i) {
            length(strsplit(i, '[.]')[[1]])
        })
        n <- names(which(z==min(z)))
        return(sample(n,1))
    } else {
        return(NA)
    }
})
table(final)

## color by final
table(final)
final[final %in% names(which(table(final) <= 50))] <- NA
final <- factor(final)
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=final,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1)









############# Just plot excitatory and inhibitory

ex.cells <- names(THS.annot)[grepl('Ex', THS.annot)]
ex.cells.refine <- final[ex.cells]
ex.cells.refine[ex.cells.refine %in% names(which(table(ex.cells.refine) < 300))] <- NA
table(ex.cells.refine)

THS.r$plotEmbedding(type='PCA',show.legend=F,groups=ex.cells.refine,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

in.cells <- names(THS.annot)[grepl('In', THS.annot)]
in.cells.refine <- final[in.cells]
in.cells.refine[in.cells.refine %in% names(which(table(in.cells.refine) < 300))] <- NA
table(in.cells.refine)

THS.r$plotEmbedding(type='PCA',show.legend=F,groups=in.cells.refine,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

THS.annot.refine <- as.character(THS.annot)
names(THS.annot.refine) <- names(THS.annot)
THS.annot.refine[names(ex.cells.refine)] <- ex.cells.refine
THS.annot.refine[names(in.cells.refine)] <- in.cells.refine
THS.r$plotEmbedding(type='PCA',show.legend=F,groups=THS.annot.refine,mark.clusters=T,mark.cluster.cex=1, main='scTHS-seq',embeddingType='tSNE', alpha=0.1, shuffle.colors=TRUE)

table(THS.annot.refine)
refine.dg <- THS.r$getDifferentialGenes(upregulated.only = TRUE, groups=THS.annot.refine, z.threshold=5)
lvec <- do.call(rbind, lapply(refine.dg, function(g) {
    g <- names(g)
    lv <- Matrix::rowSums(THS.r$counts[,g])/length(g)/Matrix::rowSums(THS.r$counts)
    lv[is.na(THS.annot)] <- 0
    lv
}))
lvec.scale <- t(lvec)

par(mfrow=c(5,2))
lapply(seq_len(length(refine.dg)), function(i) {
THS.r$plotEmbedding(type='PCA',show.legend=F,colors=lvec.scale[,i], embeddingType='tSNE', main=colnames(lvec.scale)[i])
})
