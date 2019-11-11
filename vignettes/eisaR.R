## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----availableOnline, eval=FALSE-------------------------------------------
#  pkgs <- BiocManager::available("TxDb")

## ----annotation, message=FALSE---------------------------------------------
# load package
library(eisaR)

# get TxDb object
txdbFile <- system.file("extdata", "hg19sub.sqlite", package = "eisaR")
txdb <- AnnotationDbi::loadDb(txdbFile)

## ----regions---------------------------------------------------------------
# extract filtered exonic and gene body regions
regS <- getRegionsFromTxDb(txdb = txdb, strandedData = TRUE)
regU <- getRegionsFromTxDb(txdb = txdb, strandedData = FALSE)

lengths(regS)
lengths(regU)

regS$exons

## ----exportregions---------------------------------------------------------
library(rtracklayer)
export(regS$exons, "hg19sub_exons_stranded.gtf")
export(regS$genebodies, "hg19sub_genebodies_stranded.gtf")

## ----extdata---------------------------------------------------------------
library(QuasR)
file.copy(system.file(package = "QuasR", "extdata"), ".", recursive = TRUE)

## ----align-----------------------------------------------------------------
sampleFile <- "extdata/samples_chip_single.txt"
genomeFile <- "extdata/hg19sub.fa"

proj <- qAlign("extdata/samples_rna_single.txt", "extdata/hg19sub.fa",
               splicedAlignment = TRUE)
alignmentStats(proj)

## ----count-----------------------------------------------------------------
cntEx <- qCount(proj, regU$exons, orientation = "any")
cntGb <- qCount(proj, regU$genebodies, orientation = "any")
cntIn <- cntGb - cntEx
head(cntEx)
head(cntIn)

## ----loadcounts------------------------------------------------------------
cntEx <- readRDS(system.file("extdata",
                             "Fig3abc_GSE33252_rawcounts_exonic.rds",
                             package = "eisaR"))
cntIn <- readRDS(system.file("extdata",
                             "Fig3abc_GSE33252_rawcounts_intronic.rds",
                             package = "eisaR"))

## ----runEISA---------------------------------------------------------------
# remove "width" column
Rex <- cntEx[, colnames(cntEx) != "width"]
Rin <- cntIn[, colnames(cntIn) != "width"]

# create condition factor (contrast will be TN - ES)
cond <- factor(c("ES", "ES", "TN", "TN"))

# run EISA
res <- runEISA(Rex, Rin, cond)

## ----compare---------------------------------------------------------------
res1 <- runEISA(Rex, Rin, cond, method = "published")
res2 <- runEISA(Rex, Rin, cond, method = "new")

# number of quantifyable genes
nrow(res1$DGEList)
nrow(res2$DGEList)

# number of genes with significant post-transcriptional regulation
sum(res1$tab.ExIn$FDR < 0.05)
sum(res2$tab.ExIn$FDR < 0.05)

# comparison of deltas
ids <- intersect(rownames(res1$DGEList), rownames(res2$DGEList))
cor(res1$contrasts[ids,"Dex"], res2$contrasts[ids,"Dex"])
cor(res1$contrasts[ids,"Din"], res2$contrasts[ids,"Din"])
cor(res1$contrasts[ids,"Dex.Din"], res2$contrasts[ids,"Dex.Din"])
plot(res1$contrasts[ids,"Dex.Din"], res2$contrasts[ids,"Dex.Din"])

## ----plotEISA--------------------------------------------------------------
plotEISA(res, contrast = "ExIn")

## ----normalization---------------------------------------------------------
# remove column "width"
Rex <- cntEx[,colnames(cntEx) != "width"]
Rin <- cntIn[,colnames(cntIn) != "width"]
Rall <- Rex + Rin
fracIn <- colSums(Rin)/colSums(Rall)
summary(fracIn)

# scale counts to the mean library size,
# separately for exons and introns
Nex <- t(t(Rex) / colSums(Rex) * mean(colSums(Rex)))
Nin <- t(t(Rin) / colSums(Rin) * mean(colSums(Rin)))

# log transform (add a pseudocount of 8)
NLex <- log2(Nex + 8)
NLin <- log2(Nin + 8)

## ----quantgenes------------------------------------------------------------
quantGenes <- rownames(Rex)[ rowMeans(NLex) > 5.0 & rowMeans(NLin) > 5.0 ]
length(quantGenes)

## ----dIdE------------------------------------------------------------------
Dex <- NLex[,c("MmTN_RNA_total_a","MmTN_RNA_total_b")] - NLex[,c("MmES_RNA_total_a","MmES_RNA_total_b")]
Din <- NLin[,c("MmTN_RNA_total_a","MmTN_RNA_total_b")] - NLin[,c("MmES_RNA_total_a","MmES_RNA_total_b")]
Dex.Din <- Dex - Din

cor(Dex[quantGenes,1], Dex[quantGenes,2])
cor(Din[quantGenes,1], Din[quantGenes,2])
cor(Dex.Din[quantGenes,1], Dex.Din[quantGenes,2])

## ----sig-------------------------------------------------------------------
# create DGEList object with exonic and intronic counts
library(edgeR)
cnt <- data.frame(Ex = Rex, In = Rin)
y <- DGEList(counts = cnt, genes = data.frame(ENTREZID = rownames(cnt)))

# select quantifyable genes and normalize
y <- y[quantGenes, ]
y <- calcNormFactors(y)

# design matrix with interaction term
region <- factor(c("ex","ex","ex","ex","in","in","in","in"), levels = c("in", "ex"))
cond <- rep(factor(c("ES","ES","TN","TN")), 2)
design <- model.matrix(~ region * cond)
rownames(design) <- colnames(cnt)
design

# estimate model parameters
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

# calculate likelihood-ratio between full and reduced models
lrt <- glmLRT(fit)

# create results table
tt <- topTags(lrt, n = nrow(y), sort.by = "none")
head(tt$table[order(tt$table$FDR, decreasing = FALSE), ])

## ----plot------------------------------------------------------------------
sig     <- tt$table$FDR < 0.05
sum(sig)
sig.dir <- sign(tt$table$logFC[sig])
cols <- ifelse(sig, ifelse(tt$table$logFC > 0, "#E41A1CFF", "#4DAF4AFF"), "#22222244")

# volcano plot
plot(tt$table$logFC, -log10(tt$table$FDR), col = cols, pch = 20,
     xlab = expression(paste("RNA change (log2 ",Delta,"exon/",Delta,"intron)")),
     ylab = "Significance (-log10 FDR)")
abline(h = -log10(0.05), lty = 2)
abline(v = 0, lty = 2)
text(x = par("usr")[1] + 3 * par("cxy")[1], y = par("usr")[4], adj = c(0,1),
     labels = sprintf("n=%d", sum(sig.dir == -1)), col = "#4DAF4AFF")
text(x = par("usr")[2] - 3 * par("cxy")[1], y = par("usr")[4], adj = c(1,1),
     labels = sprintf("n=%d", sum(sig.dir ==  1)), col = "#E41A1CFF")

# Delta I vs. Delta E
plot(rowMeans(Din)[quantGenes], rowMeans(Dex)[quantGenes], pch = 20, col = cols,
     xlab = expression(paste(Delta,"intron (log2 TN/ES)")),
     ylab = expression(paste(Delta,"exon (log2 TN/ES)")))
legend(x = "bottomright", bty = "n", pch = 20, col = c("#E41A1CFF","#4DAF4AFF"),
       legend = sprintf("%s (%d)", c("Up","Down"), c(sum(sig.dir == 1), sum(sig.dir == -1))))

## ----sessionInfo-----------------------------------------------------------
sessionInfo()

