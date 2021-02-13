#### CAGE Enhancer candidate analysis
#### Axel Thieffry - December 2020
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(TeMPO)
library(CAGEfightR)
library(tidyquant)
library(limma)
library(edgeR)
library(RColorBrewer)
library(BiocParallel)
register(MulticoreParam(workers=4))
library(patchwork)
library(WriteXLS)
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
'h' <- head
'l' <- length
remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}



# 0. READ INPUT DATA ####
# -----------------------
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
enhancers <- readRDS('data/RDS files/SE_enhancers.rds')
design <- readRDS('data/RDS files/design.rds')
contrasts <- readRDS('data/RDS files/contrasts.rds')
coef_levels <- factor(c('t10', 't3010', 't30', 'hen2', 'rrp4', 'hen2.t10', 'hen2.3010', 'hen2.t30', 'rrp4.t30'))



# 1. DESCRIPTIVE STATISTICS ####
# ------------------------------
# 1a. number of enhancer candidates
l(enhancers) # 155

# 1b. enhancer annotation
rowRanges(enhancers) %>%
  as.data.frame() %>%
  ggplot(aes(x=txType, fill=txType)) +
         geom_bar(lwd=.3, col='black') +
         geom_text(stat='count', aes(label=..count..), angle=90, size=7, hjust=1.5) +
         cowplot::theme_cowplot() + theme(aspect.ratio=2, legend.position='none') +
         scale_y_continuous(expand=c(0, 0)) +
         labs(x='', y='Number of enhancer candidates (N=155)')



# 2. DIFFERENTIAL EXPRESSION ANALYSIS ####
# ----------------------------------------
dge <- DGEList(assay(enhancers, 'counts')) %>% calcNormFactors(method='TMM')
v   <- voom(dge, design=design, plot=F)
fit <- lmFit(v, design=design)
fit <- contrasts.fit(fit, contrasts)
eb  <- eBayes(fit, robust=TRUE)
dt  <- decideTests(eb, adjust.method='BH', lfc=1, p.value=0.05)

# 2a. get all DEGs
coefs <- as.character(coef_levels)
names(coefs) <- coefs

degs <- map(coefs, topTable, fit=eb, number=Inf, adjust.method='BH', p.value=0.05, lfc=1) %>%
        lapply(function(x) rownames_to_column(x, 'Enh_ID')) %>%
        bind_rows(.id='coef') %>%
        mutate('direction'=ifelse(logFC > 0, 'up', 'down')) %>%
        as_tibble()

degs$effect <- case_when(degs$coef %in% c('t10', 't3010', 't30') ~ 'treatment',
                         degs$coef %in% c('hen2', 'rrp4') ~ 'exosome',
                         degs$coef %in% c('hen2.t10', 'hen2.3010', 'hen2.t30', 'rrp4.t30') ~ 'interaction',
                         TRUE ~ 'error')

# 2b. barplot - amount of DE enhancers by coefficient
deByCoef <- summary(dt) %>%
            as.data.frame() %>%
            spread_(key='Var1', value='Freq') %>%
            set_colnames(c('coefficient', 'down', 'nDE', 'up')) %>%
            mutate('coefficient'=factor(coefficient, levels=coef_levels))

melt(deByCoef, id.vars='coefficient', variable.name='direction', value.name='deEnhancers') %>%
  subset(direction != 'nDE') %>%
  ggplot(aes(x=coefficient, fill=direction, y=deEnhancers)) +
         geom_bar(stat='identity', position=position_dodge(), col='black') +
         geom_text(aes(label=deEnhancers), position=position_dodge(width=.9), hjust=-.2) +
         labs(title='DE CAGE enhancer-level',
              x='Coefficients', y='DE Genes (LogFC > 1, FDR < 0.05)') +
         coord_flip() +
         cowplot::theme_cowplot() + theme(aspect.ratio=1.2) +
         scale_y_continuous(expand=c(0, 0, 0.1, 0)) +
         scale_fill_brewer(palette='Set1', direction=-1)

# 2c. make dtdf
dtdf <- dt %>%
        as.data.frame() %>%
        rownames_to_column('BC_id') %>%
        as_tibble()

# 2d. enhancers up-regulated in rrp4 or hen2
dtdf %>%
  filter(hen2==1 | rrp4 == 1) %>%
  nrow()

# 2e. enhancers up-regulated in timecourse
dtdf %>%
  filter(t10 == 1 | t3010 == 1 | t30 == 1) %>%
  nrow()



# 3. HEATMAP OF ENHANCER CANDIDATES ####
# --------------------------------------
# 3a. make annotation for heatmap
annot_row_df <- rowRanges(enhancers) %>%
                as.data.frame() %>%
                select('Bidir.'=bidirectionality, 'Annotation'=txType) %>%
                rownames_to_column('BC_id') %>%
                left_join(select(dtdf, BC_id,
                                 t10, t3010, t30,
                                 hen2, rrp4,
                                 hen2.t30, rrp4.t30), by='BC_id') %>%
                column_to_rownames('BC_id')

annot_col_df <- colData(enhancers) %>%
                as.data.frame() %>%
                select(timepoint, genotype)

genotype <- c('wt'='#66C2A5', 'hen2'='#8DA0CB', 'rrp4'='#FC8D62')
timepoint <- c('0'='lightblue1', '10'='lightcoral', '30'='indianred4')
de_status_full <- c('-1'='navy', '0'='white', '1'='red')
de_status_limited_up <- c('0'='white', '1'='red')
de_status_limited_down <- c('-1'='navy', '0'='white')

anno_col_colors <- list('genotype'=genotype, 'timepoint'=timepoint,
                        't10'=de_status_limited_down, 't3010'=de_status_full, 't30'=de_status_full,
                        'hen2'=de_status_full, 'rrp4'=de_status_full,
                        'hen2.t30'=de_status_limited_up, 'rrp4.t30'=de_status_limited_up)

# 3b. heatmap only with DE enhancers
de_enh_ids <- degs$Enh_ID %>% unique()

enhancers[names(enhancers) %in% de_enh_ids] %>%
  assay('TPM') %>%
  add(0.1) %>%
  log() %>%
  pheatmap::pheatmap(show_rownames=F, scale='row', border_color=NA,
           cutree_cols=3, cutree_rows=3, cellwidth=10, cellheight=10,
           annotation_colors=anno_col_colors,
           annotation_row=annot_row_df, annotation_col=annot_col_df,
           clustering_method='complete')


# 4. EXPORT ENHANCER DATA ####
# ----------------------------
# 4a. all enhancer info
rowRanges(enhancers) %>%
  as.data.frame() %>%
  rownames_to_column('enhancer_id') %>%
  rename('chr'='seqnames', 'pooled_TPM'='score', 'midpoint'='thick.start', ) %>%
  select(-thick.end, -thick.width, -thick.names, -TPM_pooled) %>%
  WriteXLS(ExcelFileName='data/excel_outputs/enhancer_data.xlsx', SheetNames='all_data', row.names=F, AdjWidth=T, col.names=T, BoldHeaderRow=T)

# 4b. enhancer DEA results
degs %>%
  WriteXLS(ExcelFileName='data/excel_outputs/enhancer_DEA.xlsx', SheetNames='enhancer_DEA', row.names=F, AdjWidth=T, col.names=T, BoldHeaderRow=T)

# 4c. enhancer BED file for IGV
enhancers[names(enhancers) %in% unique(degs$Enh_ID)] %>%
  rowRanges() %>%
  export.bed('data/beds/DE_enhancers.bed')