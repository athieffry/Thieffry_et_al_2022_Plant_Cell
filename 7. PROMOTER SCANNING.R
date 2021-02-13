#### SCANNING PROMOTERS WITH TF MATRICES
#### Axel Thieffry - December 2020
set.seed(42, sample.kind='Rounding')
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(CAGEfightR)
library(RColorBrewer)
library(BiocParallel)
library(patchwork)
library(pheatmap)
library(JASPAR2020)
library(ggseqlogo)
library(TFBSTools)
register(MulticoreParam(workers=6))
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)
odb <- org.At.tair.db
options(scipen=10)
'%!in%' <- function(x,y)!('%in%'(x,y))
remove_out_of_bound <- function(GR) { idx = GenomicRanges:::get_out_of_bound_index(GR)
                                      if(length(idx) != 0) { GR[-idx]}
                                      else {GR} }
'select' <- dplyr::select
'rename' <- dplyr::rename
'h' <- head
'l' <- length



# 1. READ ALL INPUTS ####
# -----------------------
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
universe <- readRDS('data/RDS files/universe_geneID.rds')$geneID
AtTFDB <- readRDS('data/RDS files/AtTFDB_parsed.rds')
genelevel <- readRDS('data/RDS files/SE_genelevel.rds')
degs <- readRDS('data/RDS files/DEGs.rds')
clusters <- readRDS('data/RDS files/DEG_clusters.rds')
TCs <- readRDS('data/RDS files/SE_TCs.rds')



# 2. GET CAGE-DEFINED PROMOTERS ####
# ----------------------------------
if(FALSE) {
  # 2a. get TAIR10 TSSs
  TAIR10_TSSs <- genes(txdb) %>% promoters(upstream=0, downstream=1, use.names=T)
      # fix seqinfo
      seqlevels(TAIR10_TSSs) <- seqlevels(myseqinfo)
      seqlengths(TAIR10_TSSs) <- seqlengths(myseqinfo)
      genome(TAIR10_TSSs) <- genome(myseqinfo)

  # 2b. only keep TAIR10 TSSs for genes detected in CAGE experiment
  TAIR10_TSSs %<>% subset(gene_id %in% universe)

  # 2c. only keep CAGE TCs associated to a gene (on its sense strand)
  TCs_with_genes_gr <- rowRanges(TCs) %>% subset(!is.na(geneID))

  # 2d. find closest CAGE TSSs associated to the same gene, on the same strand
  TAIR10_TSSs_split <- split(TAIR10_TSSs, TAIR10_TSSs$gene_id)
  TCs_with_genes_split <- split(swapRanges(TCs_with_genes_gr), TCs_with_genes_gr$geneID)
      # sanity check
      all(names(TAIR10_TSSs_split) == names(TCs_with_genes_split))

    # this is very ugly...and very long (~10 min), but it works
    o <- GRanges()
    for(i in 1:length(TAIR10_TSSs_split)) {
            TAIR10_subset <- TAIR10_TSSs_split[i] %>% unlist()
            TCs_subset <- TCs_with_genes_split[i] %>% unlist()
            idx <- nearest(TAIR10_subset, TCs_subset)
            o <- c(o, TCs_subset[idx])
            }

  names(o) <- o$geneID
  CAGE_defined_TSSs <- o
  saveRDS(CAGE_defined_TSSs, 'data/RDS files/CAGE_defined_TSSs.rds')
}

CAGE_defined_TSSs <- readRDS('data/RDS files/CAGE_defined_TSSs.rds')



# 3. GET PROMOTER SEQUENCES ####
# ------------------------------
# 3a. get gene IDs for each cluster
clusters_splitted <- clusters %>% split(paste0('c', .$cluster))
    # sanity check
    sapply(clusters_splitted, nrow)

# 3b. get promoters regions (from closest CAGE-defined TSS to 500 bp upstream)
cluster_promoters_GR <- lapply(clusters_splitted, function(x) CAGE_defined_TSSs %>%
                                                              subset(names(.) %in% x$geneID) %>%
                                                              promoters(upstream=500, downstream=0, use.names=T))

bg_promoter_GR <- CAGE_defined_TSSs %>%
                  subset(names(.) %in% universe) %>%
                  promoters(upstream=500, downstream=0, use.names=T) %>%
                  remove_out_of_bound()

    # sanity checks
    sapply(cluster_promoters_GR, l)
    all(width(bg_promoter_GR)==500)
    all(unlist(sapply(cluster_promoters_GR, function(x) width(x)==500))==TRUE)

# 3c. get promoter sequences
genome <- BSgenome.Athaliana.TAIR.TAIR9
genome@seqinfo@genome[] <- "TAIR10" # little hack

cluster_promoters_seqs <- lapply(cluster_promoters_GR, function(x) getSeq(genome, x))
bg_promoter_seqs <- getSeq(genome, bg_promoter_GR)

# 3d. export as FASTA for RSAT scanning
fasta_filenames <- paste0('data/', names(cluster_promoters_seqs), '_promoters_500bp.fasta')
mapply(function(x, y) writeXStringSet(x, filepath=y, append=F, format='fasta'), cluster_promoters_seqs, fasta_filenames)
writeXStringSet(bg_promoter_seqs, filepath='data/bg_promoters_500bp.fasta', append=F, format='fasta')

#######
####### HERE: Do FIMO analysis online (http://meme-suite.org/tools/fimo) with default parameters, then download the TSV result files.
#######


# 4. HEATMAP RESULTS RSAT PLANT MATRIX SCAN ####
# ----------------------------------------------
# 4a. read, clean, and parse results from FIMO
fimo_files <- list.files('data/JASPAR/FIMO', pattern='FIMO_C', full.names=T)
fimo_results <- lapply(fimo_files, function(x) read.table(x, h=T, sep='\t', stringsAsFactors=F) %>% as_tibble() %>% select(motif_id, motif_alt_id, sequence_name))
names(fimo_results) <- paste0('C', 1:6)

# 4b. read background FIMO scans & compute background average matches
bg_fimo <- read.table('data/JASPAR/FIMO/FIMO_bg_1E-04.tsv', h=T, sep='\t', stringsAsFactors=F) %>%
           as_tibble() %>%
           select(motif_id, motif_alt_id, sequence_name)

bg_fimo %<>%
  group_by(motif_id, motif_alt_id) %>%
  summarize('bg_matches'=n()) %>%
  mutate('mean_bg_matches'=bg_matches/length(universe)) %>%
  select(-bg_matches) %>%
  ungroup()

# 4c. normalize C1 matches to background average
fimo_matrices_normalized <- lapply(fimo_results, function(x) group_by(x, motif_id, motif_alt_id, sequence_name) %>%
                                                             summarize('n_matches'=n()) %>%
                                                             ungroup() %>%
                                                             left_join(bg_fimo, by=c('motif_id', 'motif_alt_id')) %>%
                                                             mutate('norm_matches'=n_matches/mean_bg_matches) %>%
                                                             unite('TF', c('motif_alt_id', 'motif_id'), sep=' - ') %>%
                                                             select(TF, sequence_name, norm_matches) %>%
                                                             pivot_wider(names_from='TF', values_from='norm_matches') %>%
                                                             column_to_rownames('sequence_name'))

# 4d. set NAs to zero
fimo_matrices_normalized_nozero <- fimo_matrices_normalized
fimo_matrices_normalized_nozero %<>% lapply(function(x) { x[is.na(x)] <- 0 ; x })

# 4e. make annotation dataframes
anot_cols <- data.frame(row.names=colnames(fimo_matrices_normalized_nozero$C1),
                        'TF_fam'=c(rep('AP2-EREBP', 2), 'ARF', 'C2C2-GATA', rep('AP2-EREBP', 4), 'MYB', 'HOMEOBOX', 'HSF', 'AP2-EREBP'))

# 4f. plot normalized heatmap
promoters_left <- lapply(fimo_matrices_normalized_nozero, function(x) n_distinct(rownames(x)))

pheatmap(fimo_matrices_normalized_nozero$C1, annotation_col=anot_cols, clustering_method='ward.D2', show_rownames=F, cellwidth=10, scale='row', main=paste0('C1 - TFBS enrichment (N=', promoters_left$C1,')'))
pheatmap(fimo_matrices_normalized_nozero$C2, annotation_col=anot_cols, clustering_method='ward.D2', show_rownames=F, cellwidth=10, scale='row', main=paste0('C2 - TFBS enrichment (N=', promoters_left$C2,')'))
pheatmap(fimo_matrices_normalized_nozero$C3, annotation_col=anot_cols, clustering_method='ward.D2', show_rownames=F, cellwidth=10, scale='row', main=paste0('C3 - TFBS enrichment (N=', promoters_left$C3,')'))
pheatmap(fimo_matrices_normalized_nozero$C4, annotation_col=anot_cols, clustering_method='ward.D2', show_rownames=F, cellwidth=10, scale='row', main=paste0('C4 - TFBS enrichment (N=', promoters_left$C4,')'))
pheatmap(fimo_matrices_normalized_nozero$C5, annotation_col=anot_cols, clustering_method='ward.D2', show_rownames=F, cellwidth=10, scale='row', main=paste0('C5 - TFBS enrichment (N=', promoters_left$C5,')'))
pheatmap(fimo_matrices_normalized_nozero$C6, annotation_col=anot_cols, clustering_method='ward.D2', show_rownames=F, cellwidth=10, scale='row', main=paste0('C6 - TFBS enrichment (N=', promoters_left$C6,')'))



# 5. PLOT TF LOGOS ####
# ---------------------
# 5a. get all matrix IDs
matrices <- list.files('data/JASPAR/', pattern='.pfm') %>%
            str_remove('\\.pfm') %>% str_split('_', simplify=T) %>%
            as.data.frame() %>%
            set_colnames(c('name', 'matrix_ID'))

# 5b. get matrices
pfms <- lapply(matrices$matrix_ID, function(x) as.matrix(getMatrixByID(x=JASPAR2020, ID=x)))
names(pfms) <- paste0(matrices$name, ' (', matrices$matrix_ID, ')')

# c. plot
ggseqlogo(pfms, ncol=2) + theme_bw()
