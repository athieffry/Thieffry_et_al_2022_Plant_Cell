#### CAGEfightR core analysis
#### Axel Thieffry - November 2020
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(CAGEfightR)
library(RColorBrewer)
library(BiocParallel)
library(TeMPO)
register(MulticoreParam(workers=6))
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)
odb <- org.At.tair.db
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename



# 1. PREPARATION ####
# -------------------
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
colData <- readRDS('data/RDS files/colData.rds')


# 2. MAKE CAGE CTSSs ####
# -----------------------
# 2a. Create CTSS object
bw_plus   <- list.files(path='data/CAGE bws/bigwigs_raw', pattern='*.plus.chrWchr.bw', full.names=T) %>% BigWigFileList()
bw_minus  <- list.files(path='data/CAGE bws/bigwigs_raw', pattern='*.minus.chrWchr.bw', full.names=T) %>% BigWigFileList()
names(bw_plus) <- names(bw_minus) <- colData$Name

# 2b. quantify CTSS across samples
CTSSs <- quantifyCTSSs(plusStrand=bw_plus, minusStrand=bw_minus, genome=myseqinfo, design=colData)

#### remove noise and normalize in TPM
# 2c. calculate support based on readcounts
CTSSs <- calcSupport(CTSSs, inputAssay='counts', outputColumn='CTSS_support', unexpressed=0)

# 2d. remove CTSS supported by less than 3 libraries (min group size)
CTSSs <- subset(CTSSs, CTSS_support >= 3)

# 2e. normalize to TPM: add supported_TPM assay as well as a 'totalTags' column to colData
CTSSs <- calcTPM(CTSSs, inputAssay='counts', outputAssay='TPM', outputColumn='totalTags')
    # sanitycheck: total TPM is 1 million
    assay(CTSSs, 'TPM') %>% colSums() %>% as.data.frame()

# 2f. add TPM support
CTSSs <- calcSupport(CTSSs, inputAssay='TPM', outputColumn='TPM_support', unexpressed=0)

# 2g. calculate pooled TPM: as score in the rowRanges
CTSSs <- calcPooled(CTSSs, inputAssay='TPM', outputColumn='pooled_TPM')



# 3. UNIDIRECTIONAL TAG CLUSTERING ####
# -------------------------------------
# 3a. clustering
score(rowRanges(CTSSs)) <- rowRanges(CTSSs)$pooled_TPM # score column is needed for the next step, use pooled_TPM signal
TCs <- clusterUnidirectionally(CTSSs)

# 3b. quantify TCs with CTSS counts
TCs <- quantifyClusters(CTSSs, clusters=TCs, inputAssay='counts', sparse=F)

# 3c. calculate TPM from CTSS counts
TCs <- calcTPM(TCs, totalTags='totalTags', inputAssay='counts', outputAssay='TPM')

# 3d. calculate TPM support
TCs <- calcSupport(TCs, inputAssay='TPM', outputColumn='TPM_support', unexpressed=1)

# 3e. keep only TCs with a minimum of 1 TPM in at least 3 samples (min group size)
TCs <- subset(TCs, TPM_support >= 3)



# 4. BIDIRECTIONAL TAG CLUSTERING ####
# ------------------------------------
# 4a. clustering
enhancers <- clusterBidirectionally(CTSSs, window=251)

# 4b. calculate bidirectionality support (must have counts > 0 in both arms)
enhancers <- calcBidirectionality(enhancers, samples=CTSSs)

# 4c. plot enhancer bidirectionality support
mcols(enhancers) %>%
  as.data.frame() %>%
  ggplot(aes(x=as.factor(bidirectionality))) +
        geom_bar(stat='count') +
        geom_vline(xintercept=3.5, col='red', lty=2) +
        labs(x='Support (libraries with bidirectional signal)', y='Enhancers',
             title='Bidirectional support of enhancers',
             subtitle='Sliding window length: 500bp\n') +
        theme(aspect.ratio=1)

# 4d. keep only enhancers that are bidirectional in at least in 3 samples
enhancers <- subset(enhancers, bidirectionality >= 3)

# 4e. quantify enhancers (counts)
enhancers <- quantifyClusters(CTSSs, clusters=enhancers, inputAssay='counts')

# 4f. quantify enhancers (TPM)
enhancers <- calcTPM(enhancers, totalTags='totalTags', inputAssay='counts', outputAssay='TPM')

# 4g. calculate pooled TPM as a new column
enhancers <- calcPooled(enhancers, inputAssay='TPM', outputColumn='TPM_pooled')



# 5. BUILD CUSTOM TAIR10 ANNOTATION ####
# --------------------------------------
# 5a. Build sense annotation hierarchy (first to appear is higher priority)
# NB: using a GRangesList for the txModels will make assignTxType to not use any argument/option, so everything has to be build here
if(FALSE) {
  Antis_tair <- transcripts(txdb) %>% invertStrand()
  # make reverse strand (aka PROMPT regions, antisense to gene and up to 400bp upstream of TSS) for TAIR10
  Reverse_tair <- suppressWarnings( promoters(txdb, upstream=400, downstream=0) %>% trim() %>% invertStrand() )
  # build
  custom_hierarchy_tair10 <- GRangesList('promoter' = suppressWarnings( promoters(txdb, upstream=100, downstream=100) %>% trim() %>% granges() ),
                                         'fiveUTR' = fiveUTRsByTranscript(txdb) %>% unlist() %>% granges(),
                                         'threeUTR' = threeUTRsByTranscript(txdb) %>% unlist() %>% granges(),
                                         'CDS' = cds(txdb) %>% granges(),
                                         'exon' = exons(txdb) %>% granges(),
                                         'intron' = intronsByTranscript(txdb) %>% unlist() %>% granges(),
                                         'proximal' = suppressWarnings( promoters(txdb, upstream=400, downstream=0) %>% trim() %>% granges() ),
                                         'reverse' = Reverse_tair,
                                         'antisense' = Antis_tair)
  # fix its seqinfo
  seqlevelsStyle(custom_hierarchy_tair10) <- seqlevelsStyle(myseqinfo)
  seqlevels(custom_hierarchy_tair10) <- seqlevels(myseqinfo)
  seqinfo(custom_hierarchy_tair10) <- myseqinfo
  # save
  saveRDS(custom_hierarchy_tair10, file='data/RDS files/TAIR10_custom_annotation_hierarchy.rds')
}

custom_hierarchy_tair10 <- readRDS('data/RDS files/TAIR10_custom_annotation_hierarchy.rds')

# 5b. Build antisense annotation hierarchy (first to appear is higher priority)
if(FALSE) {
  antisense_tmp <- invertStrand(custom_hierarchy_tair10)
  names(antisense_tmp) <- paste0('antisense_', names(antisense_tmp))
  # merge with sense hierarchy (after, because order defines priority)
  custom_annotation_hierarchy_extended_TAIR10 <- c(custom_hierarchy_tair10, antisense_tmp)
  # remove redundant categories
  data.frame('antisense_annot'=names(custom_annotation_hierarchy_extended_TAIR10))
  custom_annotation_hierarchy_extended_TAIR10$reverse <- NULL # because it does not overlap the gene!!!!
  custom_annotation_hierarchy_extended_TAIR10$antisense <- NULL # because it will be annotated
  custom_annotation_hierarchy_extended_TAIR10$antisense_antisense <- NULL # because redundant
  custom_annotation_hierarchy_extended_TAIR10$antisense_promoter <- NULL # because it will either be 5'UTR or reverse
  custom_annotation_hierarchy_extended_TAIR10$antisense_reverse <- NULL # because reverse already captures them
  custom_annotation_hierarchy_extended_TAIR10$antisense_proximal <- NULL # because reverse already captures them
  # save
  saveRDS(custom_annotation_hierarchy_extended_TAIR10, 'data/RDS files/TAIR10_custom_annotation_hierarchy_extended_antisense.rds')
}

custom_hierarchy_tair10_extended_antisense <- readRDS('data/RDS files/TAIR10_custom_annotation_hierarchy_extended_antisense.rds')



# 6. ANNOTATE CAGE TCs & enhancers ####
# -------------------------------------
# 6a. fix seqinfo of Arabidopsis transcripts database
seqlevelsStyle(txdb) <- seqlevelsStyle(myseqinfo)
seqlevels(txdb) <- seqlevels(myseqinfo)

# 6b. assign TX to CAGE TCs: sense TCs are up to 400bp upstream of the TX start will be assigned to that TX
TCs <- assignTxID(TCs, txModels=txdb, outputColumn='txID', swap='thick', upstream=400, downstream=0)

# 6c. assign TxType based on custom hierarchical annotation
TCs <- assignTxType(TCs, txModels=custom_hierarchy_tair10, outputColumn='txType', swap='thick', noOverlap='intergenic')

# 6d. assign extended antisense TxType
TCs <- assignTxType(TCs, txModels=custom_hierarchy_tair10_extended_antisense, outputColumn='txType_extended', swap='thick', noOverlap='intergenic')

# 6e. relevel TxType to somewhat follow a canonical gene structure
rowRanges(TCs)$txType %<>% factor(levels=c('reverse', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR', 'antisense', 'intergenic'))
rowRanges(TCs)$txType_extended %<>% factor(levels=c('reverse', 'proximal', 'promoter', 'fiveUTR', 'antisense_fiveUTR', 'CDS', 'antisense_CDS', 'exon', 'antisense_exon', 'intron', 'antisense_intron','threeUTR', 'antisense_threeUTR', 'antisense', 'intergenic'))

# 6f. annotate enhancers
enhancers <- assignTxType(enhancers, txModels=txdb, outputColumn='txType')
rowRanges(enhancers)$txType %<>% factor(levels=c('reverse', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR', 'antisense', 'intergenic'))
enhancers <- assignGeneID(enhancers, geneModels=txdb, outputColumn='geneID', upstream=0, downstream=0, swap='thick')


# 6g. only keep intronic and intergenic enhancers
enhancers %<>% subset(txType %in% c('intron', 'intergenic')) # 170 enhancers

# 6h. only keep enhancers on canonical chromosomes
seqlevels(enhancers, pruning.mode='coarse') <- paste0('Chr', 1:5) # 155

# 6i. save
saveRDS(TCs, 'data/RDS files/SE_TCs.rds')
saveRDS(enhancers, 'data/RDS files/SE_enhancers.rds')



# 7. CAGE TC SHAPE STATISTICS ####
# --------------------------------
# 7a. temporarily assign the pooled_TPM value to the score column as it will be needed for the next step:
score(rowRanges(CTSSs)) <- rowRanges(CTSSs)$pooled_TPM

TCs <- calcShape(TCs, pooled=CTSSs, outputColumn='shape_IQR', shapeFunction=shapeIQR, lower=0.1, upper=0.9)
TCs <- calcShape(TCs, pooled=CTSSs, outputColumn='shape_entropy', shapeFunction=shapeEntropy)
TCs <- calcShape(TCs, pooled=CTSSs, outputColumn='shape_mean', shapeFunction=shapeMean)

# 7b. remove tmp column
score(rowRanges(CTSSs)) <- NULL



# 8. ANNOTATION AT GENE-LEVEL ####
# --------------------------------
# 8a. normal sense genes
TCs <- assignGeneID(TCs, geneModels=txdb, outputColumn='geneID', upstream=400, downstream=0, swap='thick')

# 8b. antisense gene (must extend to the PROMPT region before inverting strand!)
extended_genes <- genes(txdb)
promoters_genes <- promoters(extended_genes, upstream=400, downstream=0) %>% trim()
start(extended_genes) <- ifelse(strand(extended_genes)=='+', start(promoters_genes) , start(extended_genes))
end(extended_genes) <- ifelse(strand(extended_genes)=='-', end(promoters_genes), end(extended_genes))
TCs <- assignGeneID(TCs, geneModels=invertStrand(extended_genes), outputColumn='geneID_anti', upstream=0, downstream=0, swap='thick')

# 8c. gene symbols
symbols <- mapIds(odb, keys=rowRanges(TCs)$geneID, keytype='TAIR', column='SYMBOL')
symbols_anti <- mapIds(odb, keys=rowRanges(TCs)$geneID_anti, keytype='TAIR', column='SYMBOL')
rowRanges(TCs)$symbol <- as.character(symbols)
rowRanges(TCs)$symbol_anti <- as.character(symbols_anti)

# 8d. Quantify expression at Gene-level
genelevel <- quantifyGenes(TCs, genes='geneID', inputAssay='counts')
genelevel <- calcTPM(genelevel, inputAssay='counts', outputAssay='TPM')



# 9. FILTERING TCs BASED ON COMPOSITION - to use for later call of DTUs ####
# --------------------------------------------------------------------------
# 9a. remove TSSs not belonging to any gene
intragenicTSSs <- subset(TCs, !is.na(geneID))

# 9b. calculate composition: the number of samples expressing TSSs above 10% of the total gene expression
intragenicTSSs <- calcComposition(intragenicTSSs, inputAssay='counts', outputColumn='composition', unexpressed=0.1, genes='geneID')

# 9c. subset for composition in a minimum of 3 libraries (smallest sample group size)
intragenicTSSs <- subset(intragenicTSSs, composition >= 3)



# 10. SAVE EVERYTHING ####
# ------------------------
saveRDS(CTSSs, file='data/RDS files/SE_CTSSs.rds')
saveRDS(TCs, file='data/RDS files/SE_TCs.rds')
saveRDS(enhancers, file='data/RDS files/SE_enhancers.rds')
saveRDS(genelevel, file='data/RDS files/SE_genelevel.rds')
saveRDS(intragenicTSSs, file='data/RDS files/SE_intragenicTSSs.rds')

# 10b. bedfiles
rowRanges(TCs) %>% export.bed('data/beds/TCs.bed')
rowRanges(enhancers) %>% export.bed('data/beds/enhancers.bed')
