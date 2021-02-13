#### TIMECOURSE CAGE DEGs CLUSTERING
#### Axel Thieffry - November 2020
set.seed(42, sample.kind='Rounding')
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(CAGEfightR)
library(RColorBrewer)
library(BiocParallel)
library(edgeR)
library(ggrepel)
library(limma)
library(eulerr)
library(DESeq2)
library(kableExtra)
library(venn)
library(patchwork)
library(gprofiler2)
library(pheatmap)
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
genelevel <- readRDS('data/RDS files/SE_genelevel.rds')
universe <- readRDS('data/RDS files/universe_geneID.rds')$geneID
degs <- readRDS('data/RDS files/DEGs.rds')

# AtTFDB (from March 11, 2019)
if(FALSE) {
    AtTFDB <- read.table('data/TFDBs/AtTFDB/families_data.txt', sep='\t', h=F, quote='"', stringsAsFactors=F, strip.white=T) %>%
      set_colnames(c('FamilyID', 'LocusName', 'GeneName', 'Description', 'RNIntegration', 'SubFamily', 'BindingSite', 'SpecialInfo', 'bla')) %>%
      as_tibble() %>%
      select(-RNIntegration, -BindingSite:-bla, -SubFamily) %>%
      set_colnames(c('family', 'geneID', 'symbol', 'description', 'subfamily')) %>%
      mutate('geneID'=toupper(geneID),
             'family'=toupper(family),
             'symbol'=toupper(symbol)) %>%
      distinct(geneID, .keep_all=T)

    saveRDS(AtTFDB, 'data/RDS files/AtTFDB_parsed.rds')
}
AtTFDB <- readRDS('data/RDS files/AtTFDB_parsed.rds')



# 2. DEGS BARPLOT ####
# --------------------
# 2a. barplot
coef_levels <- factor(c('t10', 't3010', 't30', 'hen2', 'rrp4', 'hen2.t10', 'hen2.3010', 'hen2.t30', 'rrp4.t30'))
degs %<>% mutate('coef'=factor(coef, levels=coef_levels))

metadata(genelevel)$res %>%
  summary() %>%
  as.data.frame() %>%
  set_colnames(c('direction', 'coef', 'n')) %>%
  filter(direction!='NotSig') %>%
  mutate('direction'=ifelse(direction=='Up', 'up', 'down'),
         'coef'=factor(coef, levels=coef_levels)) %>%
  ggplot(aes(x=coef, fill=direction, y=n)) +
       geom_rect(xmin=0, xmax=3.5, ymin=-Inf, ymax=Inf, fill='blue', alpha=0.006, inherit.aes=F) +
       geom_rect(xmin=3.5, xmax=5.5, ymin=-Inf, ymax=Inf, fill='grey', alpha=0.006, inherit.aes=F) +
       geom_rect(xmin=5.5, xmax=Inf, ymin=-Inf, ymax=Inf, fill='orange', alpha=0.006, inherit.aes=F) +
       geom_bar(stat='identity', position=position_dodge(), col='black') +
       geom_text(aes(label=n), position=position_dodge(width=.9), hjust=-.2) +
       scale_y_continuous(expand=c(0.01, 0, 0.15, 0)) +
       coord_flip() + theme_bw() + theme(aspect.ratio=1) +
       labs(x='Coefficients', y='Nb. DEGs', title='Differential expression at gene-level', caption='Log(FC) >= 1, FDR <= 0.05')

# 2b. GO enrichment of 90 rrp4.t30 up-regulated genes
rrp4.t30_up_genes <- degs %>% filter(coef=='rrp4.t30' & direction=='up') %>% pull(geneID) %>% unique()
go <- gost(query=rrp4.t30_up_genes, organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe)

gores <- go$result %>%
      as_tibble() %>%
      select(p_value, term_name) %>%
      arrange(p_value)

gores %>%
  mutate('p_value'=-log(p_value, base=10)) %>%
  slice_max(order_by=p_value, n=10, with_ties=T) %>%
  ggplot(aes(x=10:1, y=p_value)) +
         geom_col(fill='red', alpha=.5) +
         geom_text(aes(label=term_name, y=.5), hjust=0, col='black') +
         coord_flip() + cowplot::theme_cowplot() +
         theme(axis.text.y=element_blank(), aspect.ratio=.75) +
         scale_y_continuous(expand=c(0, 0)) +
         labs(x='Top10 GO terms', y='-log(adj.P.value, 10)', title='GSEA: rrp4.t30 up DEGs', subtitle='Up-regulated genes in rrp4.t30 (N=90)')



# 3. TIMECOURSE HEATMAP ####
# --------------------------
# 3a. keep time course DEGS only
timecourse_geneIDs <- degs %>% filter(coef %in% c('t10', 't30')) %>% pull(geneID) %>% unique()
timecourse_degs <- degs %>% filter(geneID %in% timecourse_geneIDs)

# 3b. extract expression matrices
timecourse_degs_TPM_mat <- assay(genelevel, 'TPM') %>%
                           as.data.frame() %>%
                           filter(rownames(.) %in% timecourse_degs$geneID) %>%
                           select(matches('wt')) %>%
                           as.matrix()

timecourse_degs_TPM_mat_average <- timecourse_degs_TPM_mat %>%
                                   as.data.frame() %>%
                                   rownames_to_column('geneID') %>%
                                   melt(id.vars='geneID', variable.name='sample', value.name='TPM') %>%
                                   mutate('sample'=str_remove(sample, '_R[123]')) %>%
                                   group_by(geneID, sample) %>%
                                   summarize('meanTPM'=mean(TPM)) %>%
                                   spread(sample, meanTPM) %>%
                                   column_to_rownames('geneID') %>%
                                   as.matrix()


# 3c. make annotations and colors
heatmap_colors <- colorRampPalette(c('navy', 'white', 'Red'))

annot_df <- timecourse_degs %>%
  select(coef, geneID, direction) %>%
  spread(coef, direction) %>%
  select(geneID, t30, t3010, t10) %>%
  column_to_rownames('geneID')

annot_colors <- list('t10'=c('down'='navy', 'up'='red'),
                     't3010'=c('down'='navy', 'up'='red'),
                     't30'=c('down'='navy', 'up'='red'))

# 3d. DEG heatmap
# pseudo +0.01, complete, correlation
pheat <- timecourse_degs_TPM_mat %>%
   add(0.01) %>%
   log(base=2) %>%
   pheatmap(cluster_rows=T, cluster_cols=F, show_rownames=F,
            cellwidth=20, gaps_col=c(3, 6), scale='row',
            annotation_row=annot_df, annotation_colors=annot_colors,
            color=heatmap_colors(40), cutree_rows=6,
            clustering_method='complete', clustering_distance_rows='correlation',
            main='log2(TPM+0.01) row-scaled, complete & correlation')

# 3e. get cluster memberships and add scaled TPM matrix
orders <- data.frame('geneID'=pheat$tree_row$labels,
                     'order'=pheat$tree_row$order,
                     'cluster'=cutree(pheat$tree_row, k=6)) %>% as_tibble()

timecourse_degs_TPM_mat_scaled <- timecourse_degs_TPM_mat_average %>% t() %>% scale(center=T, scale=T) %>% t() %>%
                                  as.data.frame() %>% rownames_to_column('geneID')

# 3f. plot cluster-wise expression change (z-score): using native cuttree cluster memberships
orders %>%
  left_join(timecourse_degs_TPM_mat_scaled, by='geneID') %>%
  select(cluster:wt_30) %>%
  group_by(cluster) %>%
  summarize('0'=mean(wt_0),
            '10'=mean(wt_10),
            '30'=mean(wt_30),
            'n'=n()) %>%
  melt(id.vars=c('cluster', 'n'), variable.name='time', value.name='z') %>%
  mutate('cluster'=paste0('cluster ', cluster, ' (N=', n, ')'),
         'time'=as.numeric(as.character(time))) %>%
  select(-n) %>%
  ggplot(aes(x=time, y=z, group=cluster)) +
         geom_line() +
         geom_point() +
         facet_wrap(~cluster, ncol=1) +
         labs(y='z-score', x='time (min)', caption='cluster numbers from cuttree') +
         theme_bw() + scale_color_brewer(palette='Set1', direction=-1) +
         theme(aspect.ratio=1, legend.position='none', panel.grid.minor=element_blank())

# 3g. update cluster numbers to match heatmap from top to bottom, then re-plot z-scores
orders %<>%
  mutate('topbottom_cluster'=case_when(cluster == 1 ~ 4,
                                       cluster == 2 ~ 2,
                                       cluster == 3 ~ 6,
                                       cluster == 4 ~ 1,
                                       cluster == 5 ~ 3,
                                       cluster == 6 ~ 5))

saveRDS(orders, 'data/RDS files/clusters_with_orders.rds')

orders %>%
  left_join(timecourse_degs_TPM_mat_scaled, by='geneID') %>%
  select(topbottom_cluster:wt_30) %>%
  group_by(topbottom_cluster) %>%
  summarize('0'=mean(wt_0),
            '10'=mean(wt_10),
            '30'=mean(wt_30),
            'n'=n()) %>%
  melt(id.vars=c('topbottom_cluster', 'n'), variable.name='time', value.name='z') %>%
  mutate('topbottom_cluster'=paste0('C', topbottom_cluster, ' (N=', n, ')'),
         'time'=as.numeric(as.character(time))) %>%
  rename('cluster'='topbottom_cluster') %>%
  select(-n) %>%
  ggplot(aes(x=time, y=z, group=cluster)) +
         geom_line() +
         geom_point() +
         facet_wrap(~cluster, ncol=1) +
         labs(y='z-score', x='time (min)', caption='cluster numbers from top to bottom of heatmap') +
         theme_bw() + scale_color_brewer(palette='Set1', direction=-1) +
         theme(aspect.ratio=1, legend.position='none', panel.grid.minor=element_blank())

# 3h. GO enrichment analyses
clusters <- orders %>% select(-cluster, -order, 'cluster'=topbottom_cluster)

gobp_clusters <- split(x=clusters$geneID, f=clusters$cluster) %>%
                 lapply(function(x) gost(query=x, organism='athaliana', custom_bg=universe, sources='GO:BP',
                                         significant=T, user_threshold=0.05, correction_method='fdr'))
gomf_clusters <- split(x=clusters$geneID, f=clusters$cluster) %>%
                 lapply(function(x) gost(query=x, organism='athaliana', custom_bg=universe, sources='GO:MF',
                                         significant=T, user_threshold=0.05, correction_method='fdr'))

gobp_clusters %<>% lapply(function(x) arrange(select(as_tibble(x$result), p_value, term_name), p_value))
gomf_clusters %<>% lapply(function(x) arrange(select(as_tibble(x$result), p_value, term_name), p_value))

gobp_clusters %<>% lapply(function(x) x %>%
                                      mutate('p_value'=-log(p_value, base=10)) %>%
                                      slice_max(order_by=p_value, n=10, with_ties=F))
gomf_clusters %<>% lapply(function(x) x %>%
                                      mutate('p_value'=-log(p_value, base=10)) %>%
                                      slice_max(order_by=p_value, n=10, with_ties=F))

gobp_clusters %<>% plyr::ldply(.id='cluster')
gomf_clusters %<>% plyr::ldply(.id='cluster')

ggplot(gobp_clusters, aes(x=as.factor(60:1), y=p_value)) +
       geom_col(fill='darkgreen', alpha=.5) +
       geom_text(aes(label=term_name, y=1), hjust=0, col='black') +
       facet_grid(paste0('C', cluster)~., scales='free', drop=T) +
       coord_flip() + cowplot::theme_cowplot() +
       theme(axis.text.y=element_blank()) + theme(aspect.ratio=.3) +
       scale_y_continuous(expand=c(0, 0)) +
       labs(x='Top10 GO terms', y='', title='GSEA (GO:BP)')

ggplot(gomf_clusters, aes(x=as.factor(59:1), y=p_value)) +
       geom_col(fill='red3', alpha=.5) +
       geom_text(aes(label=term_name, y=.5), hjust=0, col='black') +
       facet_grid(paste0('C', cluster)~., scales='free', drop=T) +
       coord_flip() + cowplot::theme_cowplot() +
       theme(axis.text.y=element_blank()) + theme(aspect.ratio=.3) +
       scale_y_continuous(expand=c(0, 0)) +
       labs(x='Top10 GO terms', y='', title='GSEA (GO:MF)')

# 3i. known TFs in clusters
clusters %<>%
  mutate('isTF'=ifelse(geneID %in% AtTFDB$geneID, T, F)) %>%
  left_join(AtTFDB, by='geneID')

saveRDS(clusters, 'data/RDS files/DEG_clusters.rds')

  # absolute number of known TFs per cluster
  clusters %>%
    mutate('cluster'=factor(cluster, levels=6:1)) %>%
    ggplot(aes(x=cluster, fill=isTF)) +
           geom_bar(col='black', lwd=.3, alpha=.75) +
           geom_text(stat='count', aes(label=..count..), position=position_stack(vjust=.5)) +
           cowplot::theme_cowplot() + theme(aspect.ratio=.3) +
           scale_y_continuous(expand=c(0, 0)) +
           labs(title='TFs in clusters', caption='compared to AtTFDB', x='Clusters', y='Number of DEGs') +
           scale_fill_brewer(palette='Set2', direction=-1, name='known TF') +
           coord_flip()

  # proportion of TFs (by cluster size)
  clusters %>%
    group_by(cluster) %>%
    summarize('n_TFs'=sum(isTF),
              'n_tot'=n()) %>%
    mutate('prop'=n_TFs/n_tot,
           'cluster'=factor(cluster, levels=6:1)) %>%
    ggplot(aes(x=as.factor(cluster), y=prop, label=paste0('N=',n_TFs))) +
           geom_col() +
           geom_text(col='white', hjust=1.5) +
           cowplot::theme_cowplot() + theme(aspect.ratio=.3) +
           scale_y_continuous(expand=c(0,0)) +
           labs(title='Proportion of TFs per DEG cluster',
                x='Clusters', y='Proportion of known TFs (%)') +
           coord_flip()

  # TF number vs proportion: scatter plot
  clusters %>%
    group_by(cluster) %>%
    summarize('n_TFs'=sum(isTF),
              'n_tot'=n()) %>%
    mutate('prop'=n_TFs/n_tot,
           'cluster'=factor(cluster, levels=6:1)) %>%
    ggplot(aes(x=n_TFs, y=prop)) +
           geom_label(aes(label=cluster), fill='lightgreen', alpha=.2) +
           cowplot::theme_cowplot() + theme(aspect.ratio=1) +
           labs(x='Number of known TFs', y='Proportion of TFs to cluster size (%)')

  # top TF families in all DEGs
  top10_deg_tfam <- dplyr::count(clusters, family) %>%
                    na.omit() %>%
                    slice_max(n=10, order_by=n, with_ties=F) %>%
                    mutate('family'=factor(family, levels=rev(family)))

  ggplot(top10_deg_tfam, aes(x=family, y=n)) +
         geom_col() +
         geom_text(aes(label=n), col='white', hjust=1.5) +
         coord_flip() +
         cowplot::theme_cowplot() + theme(aspect.ratio=.7) +
         scale_y_continuous(expand=c(0,0)) +
         labs(title='Top 10 TF families in DEGs', x='TF family', y='Number of DEGs in timecourse')

  # top 10 TF families by clusters
  clusters %>%
    filter(!is.na(family)) %>%
    dplyr::count(family, by=cluster) %>%
    filter(family %in% top10_deg_tfam$family) %>%
    mutate('family'=factor(family, levels=rev(top10_deg_tfam$family))) %>%
    ggplot(aes(x=as.factor(by), y=n, fill=family)) +
           geom_col(col='black', lwd=.3, alpha=.7) +
           scale_fill_brewer(palette='Paired') +
           scale_y_continuous(expand=c(0.01,0)) +
           cowplot::theme_cowplot() + theme(aspect.ratio=3) +
           labs(x='Cluster', y='Nb. DEGs', title='Top 10 TF families')

  clusters %>%
    filter(!is.na(family)) %>%
    dplyr::count(family, by=cluster) %>%
    filter(family %in% top10_deg_tfam$family) %>%
    mutate('family'=factor(family, levels=top10_deg_tfam$family)) %>%
    ggplot(aes(x=family, y=n, fill=as.factor(by))) +
           geom_col(col='black', lwd=.3, alpha=0.75, position=position_fill()) +
           scale_fill_brewer(palette='Set3', name='Cluster') +
           scale_y_continuous(expand=c(0.01,0)) +
           cowplot::theme_cowplot() + theme(aspect.ratio=.5) +
           coord_flip() +
           labs(x='Top 10 TF families', y='Nb. DEGs')