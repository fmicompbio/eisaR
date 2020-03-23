## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----availableOnline, eval=FALSE----------------------------------------------
#  pkgs <- c(BiocManager::available("TxDb")
#            BiocManager::available("EnsDb"))

## ----annotation, message=FALSE------------------------------------------------
# load package
library(eisaR)

# get TxDb object
txdbFile <- system.file("extdata", "hg19sub.sqlite", package = "eisaR")
txdb <- AnnotationDbi::loadDb(txdbFile)

## ----regions------------------------------------------------------------------
# extract filtered exonic and gene body regions
regS <- getRegionsFromTxDb(txdb = txdb, strandedData = TRUE)
regU <- getRegionsFromTxDb(txdb = txdb, strandedData = FALSE)

lengths(regS)
lengths(regU)

regS$exons

## ----exportregions------------------------------------------------------------
library(rtracklayer)
export(regS$exons, "hg19sub_exons_stranded.gtf")
export(regS$genebodies, "hg19sub_genebodies_stranded.gtf")

## ----extdata------------------------------------------------------------------
library(QuasR)
file.copy(system.file(package = "QuasR", "extdata"), ".", recursive = TRUE)

## ----align--------------------------------------------------------------------
sampleFile <- "extdata/samples_chip_single.txt"
genomeFile <- "extdata/hg19sub.fa"

proj <- qAlign(sampleFile = "extdata/samples_rna_single.txt", 
               genome = "extdata/hg19sub.fa",
               aligner = "Rhisat2", splicedAlignment = TRUE)
alignmentStats(proj)

## ----count--------------------------------------------------------------------
cntEx <- qCount(proj, regU$exons, orientation = "any")
cntGb <- qCount(proj, regU$genebodies, orientation = "any")
cntIn <- cntGb - cntEx
head(cntEx)
head(cntIn)

## ----loadcounts---------------------------------------------------------------
cntEx <- readRDS(system.file("extdata",
                             "Fig3abc_GSE33252_rawcounts_exonic.rds",
                             package = "eisaR"))
cntIn <- readRDS(system.file("extdata",
                             "Fig3abc_GSE33252_rawcounts_intronic.rds",
                             package = "eisaR"))

## ----runEISA------------------------------------------------------------------
# remove "width" column
Rex <- cntEx[, colnames(cntEx) != "width"]
Rin <- cntIn[, colnames(cntIn) != "width"]

# create condition factor (contrast will be TN - ES)
cond <- factor(c("ES", "ES", "TN", "TN"))

# run EISA
res <- runEISA(Rex, Rin, cond)

