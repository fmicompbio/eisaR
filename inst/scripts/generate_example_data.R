## Genome
genome <- Biostrings::DNAStringSet(
    c(chr1 = "TTAACATTCGCTGGGGGAGATGACGAGACTAGCCGCCGCGTGGTCCTGCCGCATTATACGTGTTCAAGCGCCTACGTGGGTTGGGCAACCCGTGCCTATGGAGGCATGGACAAATTAGGTTCAACTTCAGCTACGTACGAGACCTAGAGGTAATAAGGGTATTTTACTCGGAGCATGTTTCAGTACGAACGTTAGATATC",
      chr2 = "CTATCGAAGTGGAATCTTGAAGAGCCCATCGGTTAAGGTCTCTCCAATGTCCAGCCTATTCTATGGCACGGCAGACCCGTTGTGCATCCACAGTGATAACTTACTTGGGCTCTTAATAGAGGAGTGTTGCCATTTTATCGGCTTGCACTCCAATTAGCACCAAGTGCCGTTATTGGGGTATTGCACTCATCAATAGCGTG")
)
genomeRevcompl <- Biostrings::reverseComplement(genome)

## Transcript-to-gene mapping
t2g <- data.frame(
    transcript_id = c("tx1.1", "tx1.2", "tx1.3", "tx2.1-I", "tx2.2-U", "tx3.1", "tx3.2", "tx3.3"),
    gene_id = c("g1", "g1", "g1", "g2-I", "g2-I", "g3-U", "g3-U", "g3-U"),
    stringsAsFactors = FALSE
)

## Transcripts
txome <- Biostrings::DNAStringSet(
    c(tx1.1 = paste0(substr(genome[["chr1"]], 116L, 120L),
                     substr(genome[["chr1"]], 141L, 160L)),
      tx1.2 = paste0(substr(genome[["chr1"]], 61L, 80L),
                     substr(genome[["chr1"]], 101L, 120L),
                     substr(genome[["chr1"]], 141L, 160L)),
      tx1.3 = paste0(substr(genome[["chr1"]], 61L, 90L)),
      `tx2.1-I` = paste0(substr(genomeRevcompl[["chr1"]], 131L, 150L),
                         substr(genomeRevcompl[["chr1"]], 171L, 190L)),
      `tx2.2-U` = paste0(substr(genomeRevcompl[["chr1"]], 121L, 140L),
                         substr(genomeRevcompl[["chr1"]], 156L, 165L),
                         substr(genomeRevcompl[["chr1"]], 181L, 190L)),
      tx3.1 = paste0(substr(genome[["chr2"]], 61L, 90L),
                     substr(genome[["chr2"]], 141L, 160L)),
      tx3.2 = paste0(substr(genome[["chr2"]], 61L, 80L),
                     substr(genome[["chr2"]], 101L, 120L),
                     substr(genome[["chr2"]], 141L, 160L)),
      tx3.3 = paste0(substr(genome[["chr2"]], 61L, 90L))
    )
)

## GTF
gtf <- GenomicRanges::GRanges(
    seqnames = rep(c("chr1", "chr2"), c(18L, 10L)),
    ranges = IRanges::IRanges(
        start = c(116L, 141L, 61L, 101L, 141L, 61L, 116L, 61L, 61L, 61L,
                  51L, 11L, 61L, 36L, 11L, 11L, 11L, 11L,
                  61L, 141L, 61L, 101L, 141L, 61L, 61L, 61L, 61L, 61L),
        end = c(120L, 160L, 80L, 120L, 160L, 90L, 160L, 160L, 90L, 160L,
                70L, 30L, 80L, 45L, 20L, 70L, 80L, 80L,
                90L, 160L, 80L, 120L, 160L, 90L, 160L, 160L, 90L, 160L)
    ),
    strand = rep(c("+", "-", "+"), c(10L, 8L, 10L)),
    type = c(rep("exon", 6L), rep("transcript", 3L), "gene", rep("exon", 5L),
             rep("transcript", 2L), "gene", rep("exon", 6L),
             rep("transcript", 3L), "gene"),
    gene_id = rep(c("g1", "g2-I", "g3-U"), c(10L, 8L, 10L)),
    transcript_id = c(rep("tx1.1", 2L), rep("tx1.2", 3L), "tx1.3", "tx1.1",
                      "tx1.2", "tx1.3", NA, rep("tx2.1-I", 2L),
                      rep("tx2.2-U", 3L), "tx2.1-I", "tx2.2-U", NA,
                      rep("tx3.1", 2L), rep("tx3.2", 3L), "tx3.3", "tx3.1",
                      "tx3.2", "tx3.3", NA),
    exon_id = c("E1", "E2", "E3", "E4", "E5", "E6", NA, NA, NA, NA,
                "E7", "E8", "E9", "E10", "E11", NA, NA, NA,
                "E12", "E13", "E14", "E15", "E16", "E17", NA, NA, NA, NA)
)

rtracklayer::export(gtf, con = "small_example.gtf", format = "gtf")
Biostrings::writeXStringSet(txome, filepath = "small_example_txome.fa")
Biostrings::writeXStringSet(genome, filepath = "small_example_genome.fa")
write.table(t2g, file = "small_example_t2g.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t")
