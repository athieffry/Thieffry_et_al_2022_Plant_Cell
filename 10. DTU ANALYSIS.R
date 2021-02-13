#### DIFFSPLICE & DTU ANALYSIS
#### Axel Thieffry - January 2021
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(CAGEfightR)
library(RColorBrewer)
library(GenomicRanges)
library(limma)
library(edgeR)
library(SummarizedExperiment)
library(eulerr)
library(BiocParallel)
library(patchwork)
library(WriteXLS)
register(MulticoreParam(workers=6))
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
'count' <- dplyr::count
'l' <- length
'h' <- head



# 1. READ DATA ####
# -----------------
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
idmapping <- readRDS('data/RDS files/id_mapping.rds')
idmapping2 <- readRDS('data/RDS files/id_mapping2.rds')
design <- readRDS('data/RDS files/design.rds')
contrasts <- readRDS('data/RDS files/contrasts.rds')
intragenic_TCs <- readRDS('data/RDS files/SE_intragenicTSSs.rds')
degs <- readRDS('data/RDS files/DEG_clusters.rds')



# 2. LIMMA DIFFSPLICE ####
# ------------------------
# 2a. run limma
dge <- DGEList(assay(intragenic_TCs, 'counts')) %>% calcNormFactors(method='TMM')
v   <- voom(dge, design=design, plot=F)
fit <- lmFit(v, design=design)
fit <- contrasts.fit(fit, contrasts)
eb  <- eBayes(fit, robust=T)
ex <- diffSplice(fit=eb, geneid=rowRanges(intragenic_TCs)$geneID, exonid=names(rowRanges(intragenic_TCs)), robust=T, verbose=T)

# 2b. retrieve DTUs (~CAGE TCs), using the t-statistic (FDR-corrected)
coefs <- colnames(ex)
names(coefs) <- coefs

all_dtus <- map(coefs, topSplice, fit=ex, test='t', FDR=0.05, number=Inf) %>%
            bind_rows(.id='coef') %>%
            as_tibble()

# 2c. add direction, txType, and rename
all_dtus %<>% mutate('direction'=ifelse(logFC < 0, 'down', 'up')) %>% rename('TC_id'='ExonID', 'geneID'='GeneID')
all_dtus %<>% left_join(rowRanges(intragenic_TCs) %>% as.data.frame() %>% select(txType) %>% rownames_to_column('TC_id'), by='TC_id')

# 2d. make timecourse-only subset of DTUs
dtus_timecourse <- all_dtus %>% filter(coef %in% c('t10', 't3010', 't30'))

# 2e. filter for effect size (at least 2-fold in 1 time course coefficient)
min_fold <- 2
min_conditions <- 1

effect_sizes_timecourse <- ex$coefficients %>% as.data.frame() %>% select(t10, t3010, t30)
effect_sizes_timecourse_filtered <- effect_sizes_timecourse[apply(effect_sizes_timecourse, 1, function(x) sum(abs(x) >= min_fold) >= min_conditions), ]
effect_sizes_timecourse_filtered_TCids <- rownames(effect_sizes_timecourse_filtered) %>% unique()

dtus_timecourse %<>% filter(TC_id %in% effect_sizes_timecourse_filtered_TCids)



# 3. DESCRIPTIVE STATS (TIMECOURSE) ####
# --------------------------------------
# 3a. total DTU TCs in timecourse?
nrow(dtus_timecourse) # 708
# 3b. number of unique DTUs?
n_distinct(dtus_timecourse$TC_id) # 432
# 3c. unique number of genes with a timecourse DTU?
n_distinct(dtus_timecourse$geneID) # 225
# 3d. annotation of DTUs?
dtus_timecourse %>%
  group_by(coef, direction, txType) %>%
  summarize('n_DTUs'=n()) %>%
  ungroup() %>%
  mutate('coef'=factor(coef, levels=c('t10', 't3010', 't30'))) %>%
  ggplot(aes(x=txType, y=n_DTUs, fill=direction)) +
         geom_col(lwd=.3, col='black', position=position_dodge()) +
         geom_text(aes(label=n_DTUs), position=position_dodge(width=.9), hjust=-.15) +
         facet_grid(coef~.) +
         scale_y_continuous(expand=c(0, 0, 0.1, 0)) +
         coord_flip() + cowplot::theme_cowplot() +
         labs(x='DTU annotation', y='Number of DTU TCs in time course', title='Annotation of CAGE DTUs in time course')
# 3e. overlap of DTUs in time course?
data.frame('TC_id'=unique(dtus_timecourse$TC_id)) %>%
  mutate('t10'=ifelse(TC_id %in% filter(dtus_timecourse, coef=='t10')$TC_id, 1, 0),
         't3010'=ifelse(TC_id %in% filter(dtus_timecourse, coef=='t3010')$TC_id, 1, 0),
         't30'=ifelse(TC_id %in% filter(dtus_timecourse, coef=='t30')$TC_id, 1, 0)) %>%
  column_to_rownames('TC_id') %>%
  eulerr::euler() %>%
  plot(quantities=T, main='DTUs TCs in time course')



# 4. OUTPUT ####
# --------------
# 4a. excel
dtus_timecourse %>%
  left_join(idmapping2, by='geneID') %>%
  select(TC_id, txType, coef, t, logFC, direction, P.Value, FDR, geneID, symbol, description) %>%
  arrange(geneID) %>%
  WriteXLS('data/excel_outputs/SUPPDATA/SUPPDATA - CAGE DTUs in timecourse.xlsx', verbose=T, row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)

# 4b. RDS
saveRDS(all_dtus, 'data/RDS files/all_DTUs.rds')
saveRDS(dtus_timecourse, 'data/RDS files/DTUs_timecourse.rds')
