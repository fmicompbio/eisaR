context("Reference generation")

library(BSgenome)

gtf <- system.file("extdata", "small_example.gtf", package = "eisaR")
txome <- system.file("extdata", "small_example_txome.fa", package = "eisaR")
genome <- system.file("extdata", "small_example_genome.fa", package = "eisaR")
t2g <- system.file("extdata", "small_example_t2g.txt", package = "eisaR")

test_that("reference extraction fails with the wrong inputs", {
    expect_error(getFeatureRanges(gtf = "nonexistent"),
                 regexp = "'gtf' must be a character scalar")
    expect_error(getFeatureRanges(gtf = 1L),
                 regexp = "'gtf' must be a character scalar")
    expect_error(getFeatureRanges(gtf = c(gtf, gtf)),
                 regexp = "'gtf' must be a character scalar")

    expect_error(getFeatureRanges(gtf = gtf, featureType = "wrongvalue"),
                 regexp = "'featureType' must be a subset of")
    expect_error(getFeatureRanges(gtf = gtf, featureType = 1L),
                 regexp = "'featureType' must be a subset of")

    ## Check intron-related arguments if 'intron' is in featureType
    expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                  intronType = "wrongvalue"),
                 regexp = "'intronType' must be either")
    expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                  intronType = 1.0),
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
    expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                  intronType = "separate", flankLength = 2L,
                                  joinOverlappingIntrons = FALSE,
                                  collapseIntronsByGene = "TRUE"),
                 regexp = "'collapseIntronsByGene' must be a logical")
    expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                  intronType = "separate", flankLength = 2L,
                                  joinOverlappingIntrons = FALSE,
                                  collapseIntronsByGene = c(TRUE, FALSE)),
                 regexp = "'collapseIntronsByGene' must be a logical")
    expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                  intronType = "separate", flankLength = 2L,
                                  joinOverlappingIntrons = FALSE,
                                  keepIntronInFeature = "TRUE"),
                 regexp = "'keepIntronInFeature' must be a logical")
    expect_error(getFeatureRanges(gtf = gtf, featureType = c("spliced", "intron"),
                                  intronType = "separate", flankLength = 2L,
                                  joinOverlappingIntrons = FALSE,
                                  keepIntronInFeature = c(TRUE, FALSE)),
                 regexp = "'keepIntronInFeature' must be a logical")

    ## If 'intron' is not in featureType, intron-related arguments are not checked
    expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                  intronType = "wrongvalue"),
                 regexp = "'verbose' must be a logical scalar")
    expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                  intronType = 1.0),
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
    expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                  intronType = "separate", flankLength = 2L,
                                  collapseIntronsByGene = "TRUE"),
                 regexp = "'verbose' must be a logical scalar")
    expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                  intronType = "separate", flankLength = 2L,
                                  collapseIntronsByGene = c(TRUE, FALSE)),
                 regexp = "'verbose' must be a logical scalar")
    expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                  intronType = "separate", flankLength = 2L,
                                  keepIntronInFeature = "TRUE"),
                 regexp = "'verbose' must be a logical scalar")
    expect_error(getFeatureRanges(gtf = gtf, featureType = "spliced", verbose = "TRUE",
                                  intronType = "separate", flankLength = 2L,
                                  keepIntronInFeature = c(TRUE, FALSE)),
                 regexp = "'verbose' must be a logical scalar")

    expect_error(getFeatureRanges(gtf = gtf, verbose = 1L),
                 regexp = "'verbose' must be a logical scalar")
    expect_error(getFeatureRanges(gtf = gtf, verbose = c(TRUE, FALSE)),
                 regexp = "'verbose' must be a logical scalar")
})

test_that("feature extraction fails if there are missing required packages", {
    # these tests assume that all non-base packages are installed in
    # .Library.site and will fail if this is not the case (e.g. on BioC builders)
    skip_on_bioc()

    # set new lib paths
    withr::with_temp_libpaths({
        # unload ns
        unloadNamespace("ensembldb")
        unloadNamespace("txdbmaker")
        unloadNamespace("GenomicFeatures")
        # test
        expect_error(getFeatureRanges(gtf = gtf))
        # reload ns
        loadNamespace("GenomicFeatures")
        loadNamespace("txdbmaker")
        loadNamespace("ensembldb")
    }, action = "prefix")

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
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
                            verbose = FALSE)
    expect_length(frs, 25L)

    expect_identical(unique(unlist(frs)$type), "exon")
    expect_named(unlist(frs), unlist(frs)$transcript_id)

    ## Spliced transcripts
    splicedIds <- metadata(frs)$corrtx$spliced
    expect_identical(sort(splicedIds), sort(tx2gene$V1))
    txgtf <- subset(rtracklayer::import(gtf), type == "exon")
    txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
    txfrs <- frs[names(txgtfl)]
    expect_named(txfrs, splicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrs), ranges(txgtfl))

    expect_identical(unique(frs$tx1.1$gene_id), "g1")
    expect_identical(unique(frs$tx1.2$gene_id), "g1")
    expect_identical(unique(frs$tx1.3$gene_id), "g1")
    expect_identical(unique(frs$`tx2.1-I`$gene_id), "g2-I")
    expect_identical(unique(frs$`tx2.2-U`$gene_id), "g2-I")
    expect_identical(unique(frs$tx3.1$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.2$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.3$gene_id), "g3-U")

    ## Unspliced transcripts
    unsplicedIds <- metadata(frs)$corrtx$unspliced
    txgtflr <- range(txgtfl)
    names(txgtflr) <- paste0(names(txgtflr), "-U")
    txfrsr <- frs[names(txgtflr)]
    expect_named(txfrsr, unsplicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrsr), ranges(txgtflr))

    expect_identical(unique(frs$`tx1.1-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.2-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.3-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx2.1-I-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx2.2-U-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx3.1-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.2-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.3-U`$gene_id), "g3-U-U")

    ## Introns
    expect_identical(ranges(frs$`tx1.1-I`), IRanges(start = 118L, end = 143L))
    expect_identical(as.character(seqnames(frs$`tx1.1-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.1-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I`), IRanges(start = 78L, end = 103L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I1`), IRanges(start = 118L, end = 143L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I1`)), "+")
    expect_identical(ranges(frs$`tx3.1-I`), IRanges(start = 88L, end = 143L))
    expect_identical(as.character(seqnames(frs$`tx3.1-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.1-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I`), IRanges(start = 78L, end = 103L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I1`), IRanges(start = 118L, end = 143L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I1`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I1`)), "+")
    expect_identical(ranges(frs$`tx2.1-I-I`), IRanges(start = 28L, end = 53L))
    expect_identical(as.character(seqnames(frs$`tx2.1-I-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.1-I-I`)), "-")
    expect_identical(ranges(frs$`tx2.2-U-I`), IRanges(start = 18L, end = 38L))
    expect_identical(as.character(seqnames(frs$`tx2.2-U-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.2-U-I`)), "-")
    expect_identical(ranges(frs$`tx2.2-U-I1`), IRanges(start = 43L, end = 63L))
    expect_identical(as.character(seqnames(frs$`tx2.2-U-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.2-U-I1`)), "-")

    expect_identical(unique(frs$`tx1.1-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx1.2-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx1.2-I1`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx2.1-I-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx2.2-U-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx2.2-U-I1`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx3.1-I`$gene_id), "g3-U-I")
    expect_identical(unique(frs$`tx3.2-I`$gene_id), "g3-U-I")
    expect_identical(unique(frs$`tx3.2-I1`$gene_id), "g3-U-I")


    ## ----------------------------------------------------------------------- ##
    ## 'Separate' intron definition, join overlapping introns
    ## ----------------------------------------------------------------------- ##
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "unspliced", "intron"),
                            intronType = "separate",
                            flankLength = 9L,
                            joinOverlappingIntrons = TRUE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
                            verbose = FALSE)
    expect_length(frs, 24L)

    expect_identical(unique(unlist(frs)$type), "exon")
    expect_named(unlist(frs), unlist(frs)$transcript_id)

    ## Spliced transcripts
    splicedIds <- metadata(frs)$corrtx$spliced
    expect_identical(sort(splicedIds), sort(tx2gene$V1))
    txgtf <- subset(rtracklayer::import(gtf), type == "exon")
    txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
    txfrs <- frs[names(txgtfl)]
    expect_named(txfrs, splicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrs), ranges(txgtfl))

    ## Unspliced transcripts
    unsplicedIds <- metadata(frs)$corrtx$unspliced
    txgtflr <- range(txgtfl)
    names(txgtflr) <- paste0(names(txgtflr), "-U")
    txfrsr <- frs[names(txgtflr)]
    expect_named(txfrsr, unsplicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrsr), ranges(txgtflr))

    ## Introns
    expect_identical(ranges(frs$`tx1.1-I`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx1.1-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.1-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I`), IRanges(start = 72L, end = 109L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I1`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I1`)), "+")
    expect_identical(ranges(frs$`tx3.1-I`), IRanges(start = 82L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx3.1-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.1-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I`), IRanges(start = 72L, end = 109L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I1`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I1`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I1`)), "+")
    expect_identical(ranges(frs$`tx2.1-I-I`), IRanges(start = 22L, end = 59L))
    expect_identical(as.character(seqnames(frs$`tx2.1-I-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.1-I-I`)), "-")
    expect_identical(ranges(frs$`tx2.2-U-I`), IRanges(start = 12L, end = 69L))
    expect_identical(as.character(seqnames(frs$`tx2.2-U-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.2-U-I`)), "-")

    ## ----------------------------------------------------------------------- ##
    ## Collapse introns by gene
    ## ----------------------------------------------------------------------- ##
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "unspliced", "intron"),
                            intronType = "separate",
                            flankLength = 3L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = TRUE,
                            keepIntronInFeature = FALSE,
                            verbose = FALSE)
    expect_length(frs, 20L)

    expect_identical(unique(unlist(frs)$type), "exon")
    expect_named(unlist(frs), unlist(frs)$transcript_id)

    ## Spliced transcripts
    splicedIds <- metadata(frs)$corrtx$spliced
    expect_identical(sort(splicedIds), sort(tx2gene$V1))
    txgtf <- subset(rtracklayer::import(gtf), type == "exon")
    txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
    txfrs <- frs[names(txgtfl)]
    expect_named(txfrs, splicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrs), ranges(txgtfl))

    expect_identical(unique(frs$tx1.1$gene_id), "g1")
    expect_identical(unique(frs$tx1.2$gene_id), "g1")
    expect_identical(unique(frs$tx1.3$gene_id), "g1")
    expect_identical(unique(frs$`tx2.1-I`$gene_id), "g2-I")
    expect_identical(unique(frs$`tx2.2-U`$gene_id), "g2-I")
    expect_identical(unique(frs$tx3.1$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.2$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.3$gene_id), "g3-U")

    ## Unspliced transcripts
    unsplicedIds <- metadata(frs)$corrtx$unspliced
    txgtflr <- range(txgtfl)
    names(txgtflr) <- paste0(names(txgtflr), "-U")
    txfrsr <- frs[names(txgtflr)]
    expect_named(txfrsr, unsplicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrsr), ranges(txgtflr))

    expect_identical(unique(frs$`tx1.1-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.2-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.3-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx2.1-I-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx2.2-U-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx3.1-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.2-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.3-U`$gene_id), "g3-U-U")

    ## Introns
    expect_identical(ranges(frs$`g1-I`), IRanges(start = 78L, end = 103L))
    expect_identical(as.character(seqnames(frs$`g1-I`)), "chr1")
    expect_identical(as.character(strand(frs$`g1-I`)), "+")
    expect_identical(ranges(frs$`g1-I1`), IRanges(start = 118L, end = 143L))
    expect_identical(as.character(seqnames(frs$`g1-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`g1-I1`)), "+")
    expect_identical(ranges(frs$`g3-U-I`), IRanges(start = 78L, end = 143L))
    expect_identical(as.character(seqnames(frs$`g3-U-I`)), "chr2")
    expect_identical(as.character(strand(frs$`g3-U-I`)), "+")
    expect_identical(ranges(frs$`g2-I-I`), IRanges(start = 18L, end = 63L))
    expect_identical(as.character(seqnames(frs$`g2-I-I`)), "chr1")
    expect_identical(as.character(strand(frs$`g2-I-I`)), "-")

    expect_identical(unique(frs$`g1-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`g1-I1`$gene_id), "g1-I")
    expect_identical(unique(frs$`g2-I-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`g3-U-I`$gene_id), "g3-U-I")

    ## ----------------------------------------------------------------------- ##
    ## Introns extending outside of feature
    ## ----------------------------------------------------------------------- ##
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "unspliced", "intron"),
                            intronType = "separate",
                            flankLength = 9L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
                            verbose = FALSE)
    expect_length(frs, 25L)

    expect_identical(unique(unlist(frs)$type), "exon")
    expect_named(unlist(frs), unlist(frs)$transcript_id)

    ## Spliced transcripts
    splicedIds <- metadata(frs)$corrtx$spliced
    expect_identical(sort(splicedIds), sort(tx2gene$V1))
    txgtf <- subset(rtracklayer::import(gtf), type == "exon")
    txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
    txfrs <- frs[names(txgtfl)]
    expect_named(txfrs, splicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrs), ranges(txgtfl))

    expect_identical(unique(frs$tx1.1$gene_id), "g1")
    expect_identical(unique(frs$tx1.2$gene_id), "g1")
    expect_identical(unique(frs$tx1.3$gene_id), "g1")
    expect_identical(unique(frs$`tx2.1-I`$gene_id), "g2-I")
    expect_identical(unique(frs$`tx2.2-U`$gene_id), "g2-I")
    expect_identical(unique(frs$tx3.1$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.2$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.3$gene_id), "g3-U")

    ## Unspliced transcripts
    unsplicedIds <- metadata(frs)$corrtx$unspliced
    txgtflr <- range(txgtfl)
    names(txgtflr) <- paste0(names(txgtflr), "-U")
    txfrsr <- frs[names(txgtflr)]
    expect_named(txfrsr, unsplicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrsr), ranges(txgtflr))

    expect_identical(unique(frs$`tx1.1-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.2-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.3-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx2.1-I-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx2.2-U-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx3.1-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.2-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.3-U`$gene_id), "g3-U-U")

    ## Introns
    expect_identical(ranges(frs$`tx1.1-I`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx1.1-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.1-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I`), IRanges(start = 72L, end = 109L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I1`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I1`)), "+")
    expect_identical(ranges(frs$`tx3.1-I`), IRanges(start = 82L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx3.1-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.1-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I`), IRanges(start = 72L, end = 109L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I1`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I1`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I1`)), "+")
    expect_identical(ranges(frs$`tx2.1-I-I`), IRanges(start = 22L, end = 59L))
    expect_identical(as.character(seqnames(frs$`tx2.1-I-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.1-I-I`)), "-")
    expect_identical(ranges(frs$`tx2.2-U-I`), IRanges(start = 12L, end = 44L))
    expect_identical(as.character(seqnames(frs$`tx2.2-U-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.2-U-I`)), "-")
    expect_identical(ranges(frs$`tx2.2-U-I1`), IRanges(start = 37L, end = 69L))
    expect_identical(as.character(seqnames(frs$`tx2.2-U-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.2-U-I1`)), "-")

    expect_identical(unique(frs$`tx1.1-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx1.2-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx1.2-I1`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx2.1-I-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx2.2-U-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx2.2-U-I1`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx3.1-I`$gene_id), "g3-U-I")
    expect_identical(unique(frs$`tx3.2-I`$gene_id), "g3-U-I")
    expect_identical(unique(frs$`tx3.2-I1`$gene_id), "g3-U-I")


    ## ----------------------------------------------------------------------- ##
    ## Keep introns within feature
    ## ----------------------------------------------------------------------- ##
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "unspliced", "intron"),
                            intronType = "separate",
                            flankLength = 9L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = TRUE,
                            verbose = FALSE)
    expect_length(frs, 25L)

    expect_identical(unique(unlist(frs)$type), "exon")
    expect_named(unlist(frs), unlist(frs)$transcript_id)

    ## Spliced transcripts
    splicedIds <- metadata(frs)$corrtx$spliced
    expect_identical(sort(splicedIds), sort(tx2gene$V1))
    txgtf <- subset(rtracklayer::import(gtf), type == "exon")
    txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
    txfrs <- frs[names(txgtfl)]
    expect_named(txfrs, splicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrs), ranges(txgtfl))

    expect_identical(unique(frs$tx1.1$gene_id), "g1")
    expect_identical(unique(frs$tx1.2$gene_id), "g1")
    expect_identical(unique(frs$tx1.3$gene_id), "g1")
    expect_identical(unique(frs$`tx2.1-I`$gene_id), "g2-I")
    expect_identical(unique(frs$`tx2.2-U`$gene_id), "g2-I")
    expect_identical(unique(frs$tx3.1$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.2$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.3$gene_id), "g3-U")

    ## Unspliced transcripts
    unsplicedIds <- metadata(frs)$corrtx$unspliced
    txgtflr <- range(txgtfl)
    names(txgtflr) <- paste0(names(txgtflr), "-U")
    txfrsr <- frs[names(txgtflr)]
    expect_named(txfrsr, unsplicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrsr), ranges(txgtflr))

    expect_identical(unique(frs$`tx1.1-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.2-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.3-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx2.1-I-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx2.2-U-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx3.1-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.2-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.3-U`$gene_id), "g3-U-U")

    ## Introns
    expect_identical(ranges(frs$`tx1.1-I`), IRanges(start = 116L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx1.1-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.1-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I`), IRanges(start = 72L, end = 109L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I`)), "+")
    expect_identical(ranges(frs$`tx1.2-I1`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx1.2-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`tx1.2-I1`)), "+")
    expect_identical(ranges(frs$`tx3.1-I`), IRanges(start = 82L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx3.1-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.1-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I`), IRanges(start = 72L, end = 109L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I`)), "+")
    expect_identical(ranges(frs$`tx3.2-I1`), IRanges(start = 112L, end = 149L))
    expect_identical(as.character(seqnames(frs$`tx3.2-I1`)), "chr2")
    expect_identical(as.character(strand(frs$`tx3.2-I1`)), "+")
    expect_identical(ranges(frs$`tx2.1-I-I`), IRanges(start = 22L, end = 59L))
    expect_identical(as.character(seqnames(frs$`tx2.1-I-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.1-I-I`)), "-")
    expect_identical(ranges(frs$`tx2.2-U-I`), IRanges(start = 12L, end = 44L))
    expect_identical(as.character(seqnames(frs$`tx2.2-U-I`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.2-U-I`)), "-")
    expect_identical(ranges(frs$`tx2.2-U-I1`), IRanges(start = 37L, end = 69L))
    expect_identical(as.character(seqnames(frs$`tx2.2-U-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`tx2.2-U-I1`)), "-")

    expect_identical(unique(frs$`tx1.1-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx1.2-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx1.2-I1`$gene_id), "g1-I")
    expect_identical(unique(frs$`tx2.1-I-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx2.2-U-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx2.2-U-I1`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`tx3.1-I`$gene_id), "g3-U-I")
    expect_identical(unique(frs$`tx3.2-I`$gene_id), "g3-U-I")
    expect_identical(unique(frs$`tx3.2-I1`$gene_id), "g3-U-I")


    ## ----------------------------------------------------------------------- ##
    ## 'Collapse' intron definition
    ## ----------------------------------------------------------------------- ##
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "unspliced", "intron"),
                            intronType = "collapse",
                            flankLength = 3L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
                            verbose = FALSE)
    expect_length(frs, 22L)

    expect_identical(unique(unlist(frs)$type), "exon")
    expect_named(unlist(frs), unlist(frs)$transcript_id)

    ## Spliced transcripts
    splicedIds <- metadata(frs)$corrtx$spliced
    expect_identical(sort(splicedIds), sort(tx2gene$V1))
    txgtf <- subset(rtracklayer::import(gtf), type == "exon")
    txgtfl <- as(split(txgtf, f = txgtf$transcript_id), "GRangesList")
    txfrs <- frs[names(txgtfl)]
    expect_named(txfrs, splicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrs), ranges(txgtfl))

    expect_identical(unique(frs$tx1.1$gene_id), "g1")
    expect_identical(unique(frs$tx1.2$gene_id), "g1")
    expect_identical(unique(frs$tx1.3$gene_id), "g1")
    expect_identical(unique(frs$`tx2.1-I`$gene_id), "g2-I")
    expect_identical(unique(frs$`tx2.2-U`$gene_id), "g2-I")
    expect_identical(unique(frs$tx3.1$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.2$gene_id), "g3-U")
    expect_identical(unique(frs$tx3.3$gene_id), "g3-U")

    ## Unspliced transcripts
    unsplicedIds <- metadata(frs)$corrtx$unspliced
    txgtflr <- range(txgtfl)
    names(txgtflr) <- paste0(names(txgtflr), "-U")
    txfrsr <- frs[names(txgtflr)]
    expect_named(txfrsr, unsplicedIds, ignore.order = TRUE)
    expect_identical(ranges(txfrsr), ranges(txgtflr))

    expect_identical(unique(frs$`tx1.1-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.2-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx1.3-U`$gene_id), "g1-U")
    expect_identical(unique(frs$`tx2.1-I-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx2.2-U-U`$gene_id), "g2-I-U")
    expect_identical(unique(frs$`tx3.1-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.2-U`$gene_id), "g3-U-U")
    expect_identical(unique(frs$`tx3.3-U`$gene_id), "g3-U-U")

    ## Introns
    expect_identical(ranges(frs$`g1-I`), IRanges(start = 88L, end = 103L))
    expect_identical(as.character(seqnames(frs$`g1-I`)), "chr1")
    expect_identical(as.character(strand(frs$`g1-I`)), "+")
    expect_identical(ranges(frs$`g1-I1`), IRanges(start = 118L, end = 143L))
    expect_identical(as.character(seqnames(frs$`g1-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`g1-I1`)), "+")
    expect_identical(ranges(frs$`g2-I-I`), IRanges(start = 28L, end = 38L))
    expect_identical(as.character(seqnames(frs$`g2-I-I`)), "chr1")
    expect_identical(as.character(strand(frs$`g2-I-I`)), "-")
    expect_identical(ranges(frs$`g2-I-I1`), IRanges(start = 43L, end = 53L))
    expect_identical(as.character(seqnames(frs$`g2-I-I1`)), "chr1")
    expect_identical(as.character(strand(frs$`g2-I-I1`)), "-")
    expect_identical(ranges(frs$`g3-U-I`), IRanges(start = 88L, end = 103L))
    expect_identical(as.character(seqnames(frs$`g3-U-I`)), "chr2")
    expect_identical(as.character(strand(frs$`g3-U-I`)), "+")
    expect_identical(ranges(frs$`g3-U-I1`), IRanges(start = 118L, end = 143L))
    expect_identical(as.character(seqnames(frs$`g3-U-I1`)), "chr2")
    expect_identical(as.character(strand(frs$`g3-U-I1`)), "+")

    expect_identical(unique(frs$`g1-I`$gene_id), "g1-I")
    expect_identical(unique(frs$`g1-I1`$gene_id), "g1-I")
    expect_identical(unique(frs$`g2-I-I`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`g2-I-I1`$gene_id), "g2-I-I")
    expect_identical(unique(frs$`g3-U-I`$gene_id), "g3-U-I")
    expect_identical(unique(frs$`g3-U-I1`$gene_id), "g3-U-I")

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
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
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
    expect_identical(as.character(txs), as.character(txseqs))

    ## Unspliced transcripts
    expect_identical(as.character(seqs[["tx1.1-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 116L, 160L)))
    expect_identical(as.character(seqs[["tx1.2-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 61L, 160L)))
    expect_identical(as.character(seqs[["tx1.3-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 61L, 90L)))
    expect_identical(as.character(seqs[["tx2.1-I-U"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 11L, 70L))))
    expect_identical(as.character(seqs[["tx2.2-U-U"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 11L, 80L))))
    expect_identical(as.character(seqs[["tx3.1-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 61L, 160L)))
    expect_identical(as.character(seqs[["tx3.2-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 61L, 160L)))
    expect_identical(as.character(seqs[["tx3.3-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 61L, 90L)))

    ## Introns
    expect_identical(as.character(seqs[["tx1.1-I"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 118L, 143L)))
    expect_identical(as.character(seqs[["tx1.2-I"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 78L, 103L)))
    expect_identical(as.character(seqs[["tx1.2-I1"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 118L, 143L)))
    expect_identical(as.character(seqs[["tx2.1-I-I"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 28L, 53L))))
    expect_identical(as.character(seqs[["tx2.2-U-I"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 18L, 38L))))
    expect_identical(as.character(seqs[["tx2.2-U-I1"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 43L, 63L))))
    expect_identical(as.character(seqs[["tx3.1-I"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 88L, 143L)))
    expect_identical(as.character(seqs[["tx3.2-I"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 78L, 103L)))
    expect_identical(as.character(seqs[["tx3.2-I1"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 118L, 143L)))

    ## ----------------------------------------------------------------------- ##
    ## 'Collapse' intron definition
    ## ----------------------------------------------------------------------- ##
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "unspliced", "intron"),
                            intronType = "collapse",
                            flankLength = 3L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
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
    expect_identical(as.character(txs), as.character(txseqs))

    ## Unspliced transcripts
    expect_identical(as.character(seqs[["tx1.1-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 116L, 160L)))
    expect_identical(as.character(seqs[["tx1.2-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 61L, 160L)))
    expect_identical(as.character(seqs[["tx1.3-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 61L, 90L)))
    expect_identical(as.character(seqs[["tx2.1-I-U"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 11L, 70L))))
    expect_identical(as.character(seqs[["tx2.2-U-U"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 11L, 80L))))
    expect_identical(as.character(seqs[["tx3.1-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 61L, 160L)))
    expect_identical(as.character(seqs[["tx3.2-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 61L, 160L)))
    expect_identical(as.character(seqs[["tx3.3-U"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 61L, 90L)))

    ## Introns
    expect_identical(as.character(seqs[["g1-I"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 88L, 103L)))
    expect_identical(as.character(seqs[["g1-I1"]]),
                     as.character(Biostrings::substr(genomeseq[["chr1"]], 118L, 143L)))
    expect_identical(as.character(seqs[["g2-I-I"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 28L, 38L))))
    expect_identical(as.character(seqs[["g2-I-I1"]]),
                     as.character(Biostrings::reverseComplement(
                         Biostrings::substr(genomeseq[["chr1"]], 43L, 53L))))
    expect_identical(as.character(seqs[["g3-U-I"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 88L, 103L)))
    expect_identical(as.character(seqs[["g3-U-I1"]]),
                     as.character(Biostrings::substr(genomeseq[["chr2"]], 118L, 143L)))
})

test_that("gtf export works", {
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "intron"),
                            intronType = "collapse",
                            flankLength = 3L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
                            verbose = FALSE)
    # arguments
    f <- file.path(tempdir(), "exported_gtf.gtf")
    expect_error(exportToGtf("wrongtype", filepath = f),
                 regexp = "'grl' must be a GRangesList")
    exportToGtf(frs, filepath = f)
    rb <- rtracklayer::import(f)

    # expected results
    expect_length(rb, 17L + 6L + 8L + 6L + 3L + 3L)
    expect_identical(ranges(subset(rb, gene_id == "g1" & type == "gene")),
                     IRanges(start = 61L, end = 160L))
    expect_identical(as.character(seqnames(subset(rb, gene_id == "g1" & type == "gene"))),
                     "chr1")
    expect_identical(as.character(strand(subset(rb, gene_id == "g1" & type == "gene"))),
                     "+")

    expect_identical(ranges(subset(rb, transcript_id == "tx1.3" & type == "transcript")),
                     IRanges(start = 61L, end = 90L))
    expect_identical(
        as.character(
            seqnames(subset(rb, transcript_id == "tx1.3" & type == "transcript"))
        ),
        "chr1")
    expect_identical(
        as.character(
            strand(subset(rb, transcript_id == "tx1.3" & type == "transcript"))
            ),
        "+")
})

test_that("gtf export fails if required packages are missing", {
    # these tests assume that all non-base packages are installed in
    # .Library.site and will fail if this is not the case (e.g. on BioC builders)
    skip_on_bioc()

    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "intron"),
                            intronType = "collapse",
                            flankLength = 3L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
                            verbose = FALSE)

    # set new lib paths
    withr::with_temp_libpaths({
        # unload ns
        unloadNamespace("ensembldb")
        unloadNamespace("txdbmaker")
        unloadNamespace("GenomicFeatures")
        unloadNamespace("BSgenome")
        unloadNamespace("rtracklayer")
        # test
        expect_error(exportToGtf(frs, tempfile()))
        # load ns
        loadNamespace("rtracklayer")
        loadNamespace("GenomicFeatures")
        loadNamespace("BSgenome")
        loadNamespace("txdbmaker")
        loadNamespace("ensembldb")
    }, action = "prefix")
})

test_that("tx2gene generation works", {
    frs <- getFeatureRanges(gtf = gtf,
                            featureType = c("spliced", "intron", "unspliced"),
                            intronType = "collapse",
                            flankLength = 3L,
                            joinOverlappingIntrons = FALSE,
                            collapseIntronsByGene = FALSE,
                            keepIntronInFeature = FALSE,
                            verbose = TRUE)
    getTx2Gene(frs, file.path(tempdir(), "tx2gene.tsv"))
    txg <- getTx2Gene(frs)

    ## Generate expected tx2gene mapping
    df <- read.delim(t2g, header = FALSE, as.is = TRUE)
    colnames(df) <- c("transcript_id", "gene_id")
    dfu <- df
    dfu$transcript_id <- paste0(dfu$transcript_id, "-U")
    dfu$gene_id <- paste0(dfu$gene_id, "-U")
    df <- rbind(df, dfu)
    df <- rbind(df, data.frame(transcript_id = c("g1-I", "g1-I1", "g2-I-I",
                                                 "g2-I-I1", "g3-U-I", "g3-U-I1"),
                               gene_id = c("g1-I", "g1-I", "g2-I-I", "g2-I-I",
                                           "g3-U-I", "g3-U-I")))
    df <- df[order(df$transcript_id, df$gene_id), ]
    rownames(df) <- NULL

    txg <- txg[order(txg$transcript_id, txg$gene_id), ]
    rownames(txg) <- NULL

    expect_identical(df, txg)
})
