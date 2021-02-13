#### RNA-seq Heatmap
#### Axel Thieffry - December 2020
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(CAGEfightR)
library(RColorBrewer)
library(BiocParallel)
library(patchwork)
library(edgeR)
library(pheatmap)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
register(MulticoreParam(workers=6))
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
'l' <- length
'h' <- head



# 1. READ DATA ####
# -----------------
clusters <- readRDS('data/RDS files/clusters_with_orders.rds') %>%
            select(-cluster) %>%
            rename('cluster'='topbottom_cluster')
gene_counts <- readRDS('data/RDS files/RNA_seq_matrix.rds')



# 2. PROCESS RNA-seq ####
# -----------------------
# 2a. get gene lengths, keep only those present in RNA-seq
gene_lengths <- genes(txdb)
gene_lengths <- tibble('geneID'=gene_lengths$gene_id,
                       'width'=width(gene_lengths)) %>%
                filter(geneID %in% rownames(gene_counts))
    # check for same order
    all(gene_lengths$geneID == rownames(gene_counts))

# 2b. normalize gene expression (CPM, logCPM, RPKM)
gene_rpkm <- rpkm(gene_counts, gene.length=gene_lengths$width) %>% as.matrix()

# 2c. sort as in CAGE clusters
clusters %<>% arrange(cluster)
gene_rpkm <- gene_rpkm[match(clusters$geneID, rownames(gene_rpkm)), ]
table(rownames(gene_rpkm) == clusters$geneID)

# 2d. make colors
hcolors <- colorRampPalette(c('navy', 'white', 'red'))(60)

# 2e. get cluster cutting places
cut_sites <- dplyr::count(clusters, cluster)$n %>% cumsum()



# 3. PLOT HEATMAPS ####
# ---------------------
gene_rpkm %>%
  as.data.frame() %>%
  add(1) %>%
  log() %>%
  pheatmap(show_rownames=F, scale='row', gaps_row=cut_sites,
           cluster_rows=F, cluster_cols=F, gaps_col=c(3, 6, 12),
           cellwidth=10, color=hcolors, main='log(RPKM+1)')

