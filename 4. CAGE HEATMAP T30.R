#### Heatmap of DE TCs at 30 minutes
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
library(edgeR)
library(limma)
library(DESeq2)
library(gprofiler2)
library(patchwork)
library(pheatmap)
register(MulticoreParam(workers=6))
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
'h' <- head
'l' <- length



# 1. READ DATA ####
# -----------------
TCs <- readRDS('data/RDS files/SE_TCs.rds')
detcs <- readRDS('data/RDS files/DE_TCs.rds')



# 2. PREPARE EXPRESSION MATRICES ####
# -----------------------------------
# 2a. get CAGE TCs that are up & down regulated at t=30 minutes
t30_up_TCs <- detcs %>% filter(coef=='t30' & direction=='up') %>% select(TC_id, logFC, adj.P.Val)
t30_down_TCs <- detcs %>% filter(coef=='t30' & direction=='down') %>% select(TC_id, logFC, adj.P.Val)

# 2b. get TPM matrix for both TC sets
sample_order <- c(paste0(rep(c('wt_', 'hen2_', 'rrp4_'), each=3), rep('0_', 9), rep(c('R1', 'R2', 'R3'), 3)),
                  paste0(rep(c('wt_', 'hen2_', 'rrp4_'), each=3), rep('30_', 9), rep(c('R1', 'R2', 'R3'), 3)))

t30_up_TPM <- assay(TCs, 'TPM') %>%
              as.data.frame() %>%
              rownames_to_column('TC_id') %>%
              filter(TC_id %in% t30_up_TCs$TC_id) %>%
              select(TC_id, all_of(sample_order)) %>%
              column_to_rownames('TC_id')

t30_down_TPM <- assay(TCs, 'TPM') %>%
                as.data.frame() %>%
                rownames_to_column('TC_id') %>%
                filter(TC_id %in% t30_down_TCs$TC_id) %>%
                select(TC_id, all_of(sample_order)) %>%
                column_to_rownames('TC_id')

# 2c. make mean TPM matrix for both TC sets
t30_up_meanTPM <- t30_up_TPM %>%
                  rownames_to_column('TC_id') %>%
                  melt(id.vars='TC_id', value.name='TPM', variable.name='sample') %>%
                  mutate('sample'=str_remove(sample, '_R[123]$')) %>%
                  group_by(TC_id, sample) %>%
                  summarise('meanTPM'=mean(TPM)) %>%
                  spread(sample, meanTPM) %>%
                  column_to_rownames('TC_id') %>%
                  select('wt_0', 'hen2_0', 'rrp4_0', 'wt_30', 'hen2_30', 'rrp4_30')

t30_down_meanTPM <- t30_down_TPM %>%
                    rownames_to_column('TC_id') %>%
                    melt(id.vars='TC_id', value.name='TPM', variable.name='sample') %>%
                    mutate('sample'=str_remove(sample, '_R[123]$')) %>%
                    group_by(TC_id, sample) %>%
                    summarise('meanTPM'=mean(TPM)) %>%
                    spread(sample, meanTPM) %>%
                    column_to_rownames('TC_id') %>%
                    select('wt_0', 'hen2_0', 'rrp4_0', 'wt_30', 'hen2_30', 'rrp4_30')



# 3. HEATMAPS ####
# ----------------
# 3a. make colors
heatmap_colors <- colorRampPalette(colors=c('navy', 'white', 'red'))

# 3b. plot using averaged TPM values
t30_up_meanTPM %>%
  add(1) %>%
  log(base=2) %>%
  pheatmap(cluster_cols=F, cluster_rows=T, show_rownames=F, show_colnames=T,
           gaps_col=3, scale='row', clustering_method='ward.D2', cellwidth=20,
           clustering_distance_rows='correlation', color=heatmap_colors(100))

t30_down_meanTPM %>%
  add(1) %>%
  log(base=2) %>%
  pheatmap(cluster_cols=F, cluster_rows=T, show_rownames=F, show_colnames=T,
           gaps_col=3, scale='row', clustering_method='ward.D2', cellwidth=20,
           clustering_distance_rows='correlation', color=heatmap_colors(100))



# 4. GO ENRICHMENT ANALYSIS ####
# ------------------------------
# 4a. correspondence between TCs and their geneID
TCs_gene_correspondences <- rowRanges(TCs) %>%
                            as.data.frame() %>%
                            select(geneID) %>%
                            filter(!is.na(geneID)) %>%
                            rownames_to_column('TC_id') %>%
                            as_tibble()

# 4b. get all genes detected in the experiment (GSEA background)
universe <- unique(TCs_gene_correspondences$geneID) # 18,934 individual genes detected

# 4c. get genes with up/down-regulated CAGE TCs
t30_up_genes <- t30_up_TCs %>%
                left_join(TCs_gene_correspondences, by='TC_id') %>%
                pull(geneID) %>%
                unique()

t30_down_genes <- t30_down_TCs %>%
                  left_join(TCs_gene_correspondences, by='TC_id') %>%
                  pull(geneID) %>%
                  unique()
# 4d. perform GO enrichment analysis
go_up <- gost(query=t30_up_genes, organism='athaliana', significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, ordered_query=F)$result %>% as_tibble()
go_down <- gost(query=t30_down_genes, organism='athaliana', significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, ordered_query=F)$result %>% as_tibble()

# 4e. plot TOP10 GO enrichment
go_up %>%
   mutate('p_value_transformed'=-log(p_value, base=10)) %>%
   slice_max(order_by=p_value_transformed, n=10, with_ties=F) %>%
   ggplot(aes(x=10:1, y=p_value_transformed)) +
          geom_col(fill='red', alpha=.5) +
          geom_text(aes(label=term_name, y=2), hjust=0, col='black') +
          coord_flip() +
          cowplot::theme_cowplot() + theme(axis.text.y=element_blank(), aspect.ratio=.75) +
          scale_y_continuous(expand=c(0, 0)) +
          labs(x='Top10 GO term', y='', title='GSEA t30 min', subtitle='Genes with up-regulated CAGE TCs') +
   go_down %>%
   mutate('p_value_transformed'=-log(p_value, base=10)) %>%
   slice_max(order_by=p_value_transformed, n=10, with_ties=F) %>%
   ggplot(aes(x=10:1, y=p_value_transformed)) +
          geom_col(fill='navy', alpha=.5) +
          geom_text(aes(label=term_name, y=.3), hjust=0, col='black') +
          coord_flip() +
          cowplot::theme_cowplot() + theme(axis.text.y=element_blank(), aspect.ratio=.75) +
          scale_y_continuous(expand=c(0, 0)) +
          labs(x='Top10 GO term', y='-log(adj.P.Value, base=10)', subtitle='Genes with down-regulated CAGE TCs') +
 plot_layout(ncol=1)



# 5. SCATTER PLOTS ####
# ---------------------
# 5a. get TPM matrix for samples at 0 and 30 minutes for up-regulated TCs at 30 min
meanTPM_scatter <- t30_up_meanTPM %>%
                   filter(rownames(.) %in% t30_up_TCs$TC_id) %>%
                   rownames_to_column('TCid')

# 5b. compute fold-change
pseudocount <- 1

scatter_df <- meanTPM_scatter %>%
              column_to_rownames('TCid') %>%
              add(pseudocount) %>%
              mutate('wt_fc'=wt_30/wt_0,
                     'hen2_fc'=hen2_30/hen2_0,
                     'rrp4_fc'=rrp4_30/rrp4_0) %>%
              select(matches('fc'))

n_TCs <- nrow(scatter_df)

ggplot(scatter_df, aes(x=wt_fc, y=hen2_fc)) +
       geom_point(alpha=.5) +
       geom_smooth(method='lm', se=F) +
       geom_abline(intercept=0, slope=1, col='red', lty=2, lwd=1) +
       scale_x_log10() +
       scale_y_log10() +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       labs(x='wt FC (wt30/wt0)', y='hen2 FC (hen2 30/0)',
            title='FC scatterplot for CAGE TCs up-regulated at t30',
            subtitle=paste0('N(TCs)=', n_TCs,' from left heatmap of Figure 1D')) +
ggplot(scatter_df, aes(x=wt_fc, y=rrp4_fc)) +
       geom_point(alpha=.5) +
       geom_smooth(method='lm', se=F) +
       geom_abline(intercept=0, slope=1, col='red', lty=2, lwd=1) +
       scale_x_log10() +
       scale_y_log10() +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       labs(x='wt FC (wt30/wt0)', y='rrp4 FC (rrp4 30/0)') +
plot_layout(ncol=1)
