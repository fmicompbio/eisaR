context("Reference generation")

library(BSgenome)

gtf <- system.file("extdata/small_example.gtf", package = "eisaR")
txome <- system.file("extdata/small_example_txome.fa", package = "eisaR")
genome <- system.file("extdata/small_example_genome.fa", package = "eisaR")
t2g <- system.file("extdata/small_example_t2g.txt", package = "eisaR")

test_that("reference extraction fails with the wrong inputs", {
  expect_error(getFeatureRanges(gtf = "nonexistent"), 
               regexp = "'gtf' must be a character scalar")
  expect_error(getFeatureRanges(gtf = 1), 
               regexp = "'gtf' must be a character scalar")
  expect_error(getFeatureRanges(gtf = c(gtf, gtf)), 
               regexp = "'gtf' must be a character scalar")
  
  expect_error(getFeatureRanges(gtf = gtf, featureType = "wrongvalue"), 
               regexp = "'featureType' must be a subset of")
  expect_error(getFeatureRanges(gtf = gtf, featureType = 1), 
               regexp = "'featureType' must be a subset of")
  
  ## Check intron-related arguments if 'intron' is in featureType
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = "wrongvalue"), 
               regexp = "'intronType' must be either")
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = 1), 
               regexp = "'intronType' must be either")
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = c("separate", "collapse")), 
               regexp = "'intronType' must be either")
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = "separate", flankLength = "string"), 
               regexp = "'flankLength' must be a numeric")
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = "separate", flankLength = -1L), 
               regexp = "'flankLength' must be a numeric")
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = "separate", flankLength = c(1L, 2L)), 
               regexp = "'flankLength' must be a numeric")
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = "separate", flankLength = 2L,
                                joinOverlappingIntrons = "TRUE"), 
               regexp = "'joinOverlappingIntrons' must be a logical")
  expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                intronType = "separate", flankLength = 2L,
                                joinOverlappingIntrons = c(TRUE, FALSE)), 
               regexp = "'joinOverlappingIntrons' must be a logical")
  
  ## If 'intron' is not in featureType, intron-related arguments are not checked
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = "wrongvalue"), 
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = 1), 
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = c("separate", "collapse")), 
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = "separate", flankLength = "string"), 
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = "separate", flankLength = -1L), 
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = "separate", flankLength = c(1L, 2L)), 
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = "separate", flankLength = 2L,
                                joinOverlappingIntrons = "TRUE"), 
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                intronType = "separate", flankLength = 2L,
                                joinOverlappingIntrons = c(TRUE, FALSE)), 
               regexp = "'verbose' must be a logical scalar")
  
  expect_error(getFeatureRanges(gtf = gtf, verbose = 1),
               regexp = "'verbose' must be a logical scalar")
  expect_error(getFeatureRanges(gtf = gtf, verbose = c(TRUE, FALSE)),
               regexp = "'verbose' must be a logical scalar")
})

test_that("feature range extraction works", {
  tx2gene <- read.delim(t2g, header = FALSE, as.is = TRUE)

  ## ----------------------------------------------------------------------- ##
  ## 'Separate' intron definition
  ## ----------------------------------------------------------------------- ##
  frs <- getFeatureRanges(gtf = gtf, 
                          featureType = c("spliced", "unspliced", "intron"), 
                          intronType = "separate", 
                          flankLength = 3L, 
                          joinOverlappingIntrons = FALSE, 
                          verbose = FALSE)
  expect_length(frs, 25L)
  
  ## Spliced transcripts
  spliced_ids <- metadata(frs)$corrtx$spliced
  expect_equal(sort(spliced_ids), sort(tx2gene$V1))
  txgtf <- subset(rtracklayer::import(gtf), type == "exon")
  txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
  txfrs <- frs[names(txgtfl)]
  expect_named(txfrs, spliced_ids, ignore.order = TRUE)
  expect_equal(ranges(txfrs), ranges(txgtfl))
  
  ## Unspliced transcripts
  unspliced_ids <- metadata(frs)$corrtx$unspliced
  txgtflr <- range(txgtfl)
  names(txgtflr) <- paste0(names(txgtflr), "-U")
  txfrsr <- frs[names(txgtflr)]
  expect_named(txfrsr, unspliced_ids, ignore.order = TRUE)
  expect_equal(ranges(txfrsr), ranges(txgtflr))
  
  ## Introns
  expect_equal(ranges(frs$`tx1.1-I`), IRanges(start = 88, end = 143))
  expect_equal(as.character(seqnames(frs$`tx1.1-I`)), "chr1")
  expect_equal(as.character(strand(frs$`tx1.1-I`)), "+")
  expect_equal(ranges(frs$`tx1.2-I`), IRanges(start = 78, end = 103))
  expect_equal(as.character(seqnames(frs$`tx1.2-I`)), "chr1")
  expect_equal(as.character(strand(frs$`tx1.2-I`)), "+")
  expect_equal(ranges(frs$`tx1.2-I1`), IRanges(start = 118, end = 143))
  expect_equal(as.character(seqnames(frs$`tx1.2-I1`)), "chr1")
  expect_equal(as.character(strand(frs$`tx1.2-I1`)), "+")
  expect_equal(ranges(frs$`tx3.1-I`), IRanges(start = 88, end = 143))
  expect_equal(as.character(seqnames(frs$`tx3.1-I`)), "chr2")
  expect_equal(as.character(strand(frs$`tx3.1-I`)), "+")
  expect_equal(ranges(frs$`tx3.2-I`), IRanges(start = 78, end = 103))
  expect_equal(as.character(seqnames(frs$`tx3.2-I`)), "chr2")
  expect_equal(as.character(strand(frs$`tx3.2-I`)), "+")
  expect_equal(ranges(frs$`tx3.2-I1`), IRanges(start = 118, end = 143))
  expect_equal(as.character(seqnames(frs$`tx3.2-I1`)), "chr2")
  expect_equal(as.character(strand(frs$`tx3.2-I1`)), "+")
  expect_equal(ranges(frs$`tx2.1-I-I`), IRanges(start = 28, end = 53))
  expect_equal(as.character(seqnames(frs$`tx2.1-I-I`)), "chr1")
  expect_equal(as.character(strand(frs$`tx2.1-I-I`)), "-")
  expect_equal(ranges(frs$`tx2.2-U-I`), IRanges(start = 18, end = 38))
  expect_equal(as.character(seqnames(frs$`tx2.2-U-I`)), "chr1")
  expect_equal(as.character(strand(frs$`tx2.2-U-I`)), "-")
  expect_equal(ranges(frs$`tx2.2-U-I1`), IRanges(start = 43, end = 63))
  expect_equal(as.character(seqnames(frs$`tx2.2-U-I1`)), "chr1")
  expect_equal(as.character(strand(frs$`tx2.2-U-I1`)), "-")
  
  ## ----------------------------------------------------------------------- ##
  ## 'Collapse' intron definition
  ## ----------------------------------------------------------------------- ##
  frs <- getFeatureRanges(gtf = gtf, 
                          featureType = c("spliced", "unspliced", "intron"), 
                          intronType = "collapse", 
                          flankLength = 3L, 
                          joinOverlappingIntrons = FALSE, 
                          verbose = FALSE)
  expect_length(frs, 22L)
  
  ## Spliced transcripts
  spliced_ids <- metadata(frs)$corrtx$spliced
  expect_equal(sort(spliced_ids), sort(tx2gene$V1))
  txgtf <- subset(rtracklayer::import(gtf), type == "exon")
  txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
  txfrs <- frs[names(txgtfl)]
  expect_named(txfrs, spliced_ids, ignore.order = TRUE)
  expect_equal(ranges(txfrs), ranges(txgtfl))
  
  ## Unspliced transcripts
  unspliced_ids <- metadata(frs)$corrtx$unspliced
  txgtflr <- range(txgtfl)
  names(txgtflr) <- paste0(names(txgtflr), "-U")
  txfrsr <- frs[names(txgtflr)]
  expect_named(txfrsr, unspliced_ids, ignore.order = TRUE)
  expect_equal(ranges(txfrsr), ranges(txgtflr))
  
  ## Introns
  expect_equal(ranges(frs$`g1-I`), IRanges(start = 88, end = 103))
  expect_equal(as.character(seqnames(frs$`g1-I`)), "chr1")
  expect_equal(as.character(strand(frs$`g1-I`)), "+")
  expect_equal(ranges(frs$`g1-I1`), IRanges(start = 118, end = 143))
  expect_equal(as.character(seqnames(frs$`g1-I1`)), "chr1")
  expect_equal(as.character(strand(frs$`g1-I1`)), "+")
  expect_equal(ranges(frs$`g2-I-I`), IRanges(start = 28, end = 38))
  expect_equal(as.character(seqnames(frs$`g2-I-I`)), "chr1")
  expect_equal(as.character(strand(frs$`g2-I-I`)), "-")
  expect_equal(ranges(frs$`g2-I-I1`), IRanges(start = 43, end = 53))
  expect_equal(as.character(seqnames(frs$`g2-I-I1`)), "chr1")
  expect_equal(as.character(strand(frs$`g2-I-I1`)), "-")
  expect_equal(ranges(frs$`g3-U-I`), IRanges(start = 88, end = 103))
  expect_equal(as.character(seqnames(frs$`g3-U-I`)), "chr2")
  expect_equal(as.character(strand(frs$`g3-U-I`)), "+")
  expect_equal(ranges(frs$`g3-U-I1`), IRanges(start = 118, end = 143))
  expect_equal(as.character(seqnames(frs$`g3-U-I1`)), "chr2")
  expect_equal(as.character(strand(frs$`g3-U-I1`)), "+")
  
})
  
test_that("feature sequence extraction works", {
  ## ----------------------------------------------------------------------- ##
  ## 'Separate' intron definition
  ## ----------------------------------------------------------------------- ##
  frs <- getFeatureRanges(gtf = gtf, 
                          featureType = c("spliced", "unspliced", "intron"), 
                          intronType = "separate", 
                          flankLength = 3L, 
                          joinOverlappingIntrons = FALSE, 
                          verbose = FALSE)
  expect_length(frs, 25L)
  expect_length(metadata(frs)$featurelist, 3L)
  expect_length(metadata(frs)$featurelist$spliced, 8L)
  expect_length(metadata(frs)$featurelist$unspliced, 8L)
  expect_length(metadata(frs)$featurelist$intron, 9L)
  
  genomeseq <- Biostrings::readDNAStringSet(genome)
  seqs <- GenomicFeatures::extractTranscriptSeqs(x = genomeseq, transcripts = frs)
  
  txs <- Biostrings::readDNAStringSet(txome)
  
  ## Spliced transcripts
  txseqs <- seqs[names(txs)]
  expect_equivalent(txs, txseqs)
  expect_equal(as.character(txs), as.character(txseqs))
  
  ## Unspliced transcripts
  expect_equal(as.character(seqs[["tx1.1-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 61, 160)))
  expect_equal(as.character(seqs[["tx1.2-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 61, 160)))
  expect_equal(as.character(seqs[["tx1.3-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 61, 90)))
  expect_equal(as.character(seqs[["tx2.1-I-U"]]),
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 11, 70))))
  expect_equal(as.character(seqs[["tx2.2-U-U"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 11, 80))))
  expect_equal(as.character(seqs[["tx3.1-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 61, 160)))
  expect_equal(as.character(seqs[["tx3.2-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 61, 160)))
  expect_equal(as.character(seqs[["tx3.3-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 61, 90)))
  
  ## Introns
  expect_equal(as.character(seqs[["tx1.1-I"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 88, 143)))
  expect_equal(as.character(seqs[["tx1.2-I"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 78, 103)))
  expect_equal(as.character(seqs[["tx1.2-I1"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 118, 143)))
  expect_equal(as.character(seqs[["tx2.1-I-I"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 28, 53))))
  expect_equal(as.character(seqs[["tx2.2-U-I"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 18, 38))))
  expect_equal(as.character(seqs[["tx2.2-U-I1"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 43, 63))))
  expect_equal(as.character(seqs[["tx3.1-I"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 88, 143)))
  expect_equal(as.character(seqs[["tx3.2-I"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 78, 103)))
  expect_equal(as.character(seqs[["tx3.2-I1"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 118, 143)))
  
  ## ----------------------------------------------------------------------- ##
  ## 'Collapse' intron definition
  ## ----------------------------------------------------------------------- ##
  frs <- getFeatureRanges(gtf = gtf, 
                          featureType = c("spliced", "unspliced", "intron"), 
                          intronType = "collapse", 
                          flankLength = 3L, 
                          joinOverlappingIntrons = FALSE, 
                          verbose = FALSE)
  expect_length(frs, 22L)
  expect_length(metadata(frs)$featurelist, 3L)
  expect_length(metadata(frs)$featurelist$spliced, 8L)
  expect_length(metadata(frs)$featurelist$unspliced, 8L)
  expect_length(metadata(frs)$featurelist$intron, 6L)
  
  genomeseq <- Biostrings::readDNAStringSet(genome)
  seqs <- GenomicFeatures::extractTranscriptSeqs(x = genomeseq, transcripts = frs)
  
  txs <- Biostrings::readDNAStringSet(txome)
  
  ## Spliced transcripts
  txseqs <- seqs[names(txs)]
  expect_equivalent(txs, txseqs)
  expect_equal(as.character(txs), as.character(txseqs))
  
  ## Unspliced transcripts
  expect_equal(as.character(seqs[["tx1.1-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 61, 160)))
  expect_equal(as.character(seqs[["tx1.2-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 61, 160)))
  expect_equal(as.character(seqs[["tx1.3-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 61, 90)))
  expect_equal(as.character(seqs[["tx2.1-I-U"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 11, 70))))
  expect_equal(as.character(seqs[["tx2.2-U-U"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 11, 80))))
  expect_equal(as.character(seqs[["tx3.1-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 61, 160)))
  expect_equal(as.character(seqs[["tx3.2-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 61, 160)))
  expect_equal(as.character(seqs[["tx3.3-U"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 61, 90)))
  
  ## Introns
  expect_equal(as.character(seqs[["g1-I"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 88, 103)))
  expect_equal(as.character(seqs[["g1-I1"]]), 
               as.character(Biostrings::substr(genomeseq[["chr1"]], 118, 143)))
  expect_equal(as.character(seqs[["g2-I-I"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 28, 38))))
  expect_equal(as.character(seqs[["g2-I-I1"]]), 
               as.character(Biostrings::reverseComplement(Biostrings::substr(genomeseq[["chr1"]], 43, 53))))
  expect_equal(as.character(seqs[["g3-U-I"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 88, 103)))
  expect_equal(as.character(seqs[["g3-U-I1"]]), 
               as.character(Biostrings::substr(genomeseq[["chr2"]], 118, 143)))
})

test_that("gtf export works", {
  frs <- getFeatureRanges(gtf = gtf, 
                          featureType = c("spliced", "intron"), 
                          intronType = "collapse", 
                          flankLength = 3L, 
                          joinOverlappingIntrons = FALSE, 
                          verbose = FALSE)
  f <- file.path(tempdir(), "exported_gtf.gtf")
  exportToGtf(frs, filepath = f)
  rb <- rtracklayer::import(f)
  
  expect_length(rb, 17L + 6L + 8L + 6L + 3L + 3L)
  expect_equal(ranges(subset(rb, gene_id == "g1" & type == "gene")), 
               IRanges(start = 61, end = 160))
  expect_equal(as.character(seqnames(subset(rb, gene_id == "g1" & type == "gene"))), 
               "chr1")
  expect_equal(as.character(strand(subset(rb, gene_id == "g1" & type == "gene"))), 
               "+")
  
  expect_equal(ranges(subset(rb, transcript_id == "tx1.3" & type == "transcript")), 
               IRanges(start = 61, end = 90))
  expect_equal(as.character(seqnames(subset(rb, transcript_id == "tx1.3" & type == "transcript"))), 
               "chr1")
  expect_equal(as.character(strand(subset(rb, transcript_id == "tx1.3" & type == "transcript"))), 
               "+")
  
})


