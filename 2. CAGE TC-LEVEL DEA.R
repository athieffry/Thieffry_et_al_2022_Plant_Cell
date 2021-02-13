#### Differential Expression Analyses (CAGE TC-level)
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
library(ggridges)
library(edgeR)
library(limma)
library(DESeq2)
register(MulticoreParam(workers=6))
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)
odb <- org.At.tair.db
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
# remove out of bound from GRanges function
remove_out_of_bound <- function(GR) {
                                    idx = GenomicRanges:::get_out_of_bound_index(GR)
                                    if(length(idx) != 0) { GR[-idx]}
                                    else {GR}
                                  }
'select' <- dplyr::select
'rename' <- dplyr::rename



# 1. READ ALL INPUTS ####
# -----------------------
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
colData <- readRDS('data/RDS files/colData.rds')
TCs <- readRDS('data/RDS files/SE_TCs.rds')
enhancers <- readRDS('data/RDS files/SE_enhancers.rds')
genelevel <- readRDS('data/RDS files/SE_genelevel.rds')



# 2. DESIGN & CONTRAST MATRICES ####
# ----------------------------------
# 2a. design matrix from colData
if(FALSE) {
    # check that wt & 0 are reference levels
    levels(colData$genotype)
    levels(colData$timepoint)
    # make design matrix
    design <- model.matrix(~ genotype * timepoint, data=colData) %>% as.data.frame()
    colnames(design) <- make.names(colnames(design)) %>% str_remove('genotype') %>% str_remove('imepoint')
    design$rrp4.t10 <- NULL # remove empty level (we don't have that sample)
    # visualize (sanity check)
    annot_design <- data.frame(row.names = rownames(design),
                               'genotype' = c(rep('hen2',9), rep('rrp4',6), rep('wt',9)) %>% as.factor(),
                               'timepoint' = c(rep(c(0, 10, 30), each=3), rep(c(0, 30), each=3), rep(c(0, 10, 30), each=3)) %>% as.factor(),
                               'replicate' = rep(c(1, 2, 3), 8) %>% as.factor())

    pheatmap::pheatmap(design, cluster_rows=F, cluster_cols=F, cellwidth=10, cellheight=10, color=c('white', 'black'),
                       main='DESIGN MATRIX ( ~ genotype * timepoint)', annotation_row=annot_design, legend=F)

    # 2b. make contrast matrix
    contrasts <- diag(ncol(design))
    colnames(contrasts) <- colnames(design)
    rownames(contrasts) <- colnames(design)
    contrasts <- contrasts[ , -1]
    contrasts <- cbind(contrasts,
                       't3010' = c(0, 0, 0, -1, 1, 0, 0, 0),
                       'hen2.3010' = c(0, 0, 0, 0, 0, -1, 1, 0))

    pheatmap::pheatmap(contrasts, cluster_rows=F, cluster_cols=F, color=c('red', 'white', 'darkgreen'), cellwidth=15, cellheight=15,
                       main='CONTRAST MATRIX', legend=F, display_numbers=T, number_color='white', number_format='%s')

    # 2c. save
    saveRDS(design, file='data/RDS files/design.rds')
    saveRDS(contrasts, file='data/RDS files/contrasts.rds')
}

design <- readRDS('data/RDS files/design.rds')
contrasts <- readRDS('data/RDS files/contrasts.rds')



# 3. LIMMA DE ANALYSIS ####
# -------------------------
# 3a. run differential expression analysis
dge <- DGEList(counts=assay(TCs, 'counts'),
               samples=colData,
               lib.size=colData(TCs)$totalTags) %>%
       calcNormFactors(method='TMM')
v <- voom(dge, design=design, plot=T)
fit <- lmFit(v, design=design)
fit <- contrasts.fit(fit, contrasts)
eb <- eBayes(fit, robust=T)
res <- decideTests(eb, adjust.method='BH', lfc=1, p.value=0.05)

# 3b. get numeric matrix of normalized expression values on the log2 scale
assay(TCs, 'E') <- v$E

# 3c. multidimensional scaling plot (MDS)
mds <- limma::plotMDS(v) ; dev.off()

mds$cmdscale.out %>%
  as.data.frame() %>%
  rownames_to_column('x') %>%
  set_colnames(c('sample', 'PC1', 'PC2')) %>%
  separate('sample', c('genotype', 'timepoint', 'rep'), sep='_') %>%
  ggplot(aes(y=-PC2, x=PC1, col=timepoint, shape=genotype)) +
         geom_point(size=4) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(title='MDS plot',
              x='Dimension 1: Leading log2 fold-change',
              y='Dimension 2: Leading log2 fold-change') +
         scale_color_manual(values=c('#71A0C5', 'orange', '#BF717F')) +
         scale_shape_manual(values=c(17, 15, 16))

# 3d. extract results for all coefficients
coefs <- colnames(res)
names(coefs) <- coefs

dets <- map(coefs, topTable, fit=eb, number=Inf, adjust.method='BH', p.value=0.05, lfc=1) %>%
  lapply(function(x) rownames_to_column(x, 'TC_id')) %>%
  bind_rows(.id='coef') %>%
  mutate('direction'=ifelse(logFC > 0, 'up', 'down')) %>%
  as_tibble()

all_de <- map(coefs, topTable, fit=eb, number=Inf, adjust.method='BH') %>%
  lapply(function(x) rownames_to_column(x, 'TC_id')) %>%
  bind_rows(.id='coef') %>%
  as_tibble()

# 3e. venn diagram
coefs
par(pty='s', mar=rep(0,4))
vennDiagram(res[ , 't30'], include=c('up', 'down'), counts.col=c('red', 'blue'), names='t30', mar=rep(0, 4))
vennDiagram(res[ , c('rrp4', 'hen2')], include=c('up', 'down'), counts.col=c('red', 'blue'), names=c('rrp4', 'hen2-4'), mar=rep(0, 4))
vennDiagram(res[ , c('hen2.t30', 'rrp4.t30')], include=c('up', 'down'), counts.col=c('red', 'blue'), names=c('hen2-4 t30', 'rrp4 t30'), mar=rep(0, 4))
dev.off()



# 4. DEA POST-PROCESSING ####
# ---------------------------
# 4a. make coefficient factor level
coef_levels <- factor(c('t10', 't3010', 't30', 'hen2', 'rrp4', 'hen2.t10', 'hen2.3010', 'hen2.t30', 'rrp4.t30'))

# 4b. augment results with annotation, geneID, add effects, etc...
dets %<>% left_join(rowData(TCs) %>%
                      as.data.frame() %>%
                      rename('pooled_TPM'='score', 'CAGE_TC_peak'='thick.start', 'TC_id'='thick.names') %>%
                      select(-thick.width, -thick.end, everything()),
                    by='TC_id')

# 4c. add effect (flg22, exosome, interaction)
dets %<>% mutate('effect'=case_when(.$coef %in% c('hen2', 'rrp4') ~ 'exosome',
                                    .$coef %in% c('t10', 't30', 't3010') ~ 'flg22',
                                    .$coef %in% c('rrp4.t30') ~ 'interaction',
                                    TRUE ~ 'ERROR'))
all_de %<>% mutate('effect'=case_when(.$coef %in% c('hen2', 'rrp4') ~ 'exosome',
                                      .$coef %in% c('t10', 't30', 't3010') ~ 'flg22',
                                      .$coef %in% c('rrp4.t30', 'hen2.3010', 'hen2.t10', 'hen2.t30') ~ 'interaction',
                                      TRUE ~ 'ERROR'))

    # sanity check
    with(dets, table(coef, effect))
    with(all_de, table(coef, effect))

# 4d. add DE summary results to TCs SE
metadata(TCs)$res <- as.matrix(res)

# 4e. plot number of up/down CAGE TCs by coefficient
summary(res) %>%
  as.data.frame() %>%
  spread_(key='Var1', value='Freq') %>%
  set_colnames(c('coefficient', 'down', 'nDE', 'up')) %>%
  melt(id.vars='coefficient', variable.name='direction', value.name='deTSSs') %>%
  subset(direction != 'nDE') %>%
  mutate('coefficient'=factor(coefficient, levels=coef_levels)) %>%
  ggplot(aes(x=coefficient, fill=direction, y=deTSSs)) +
         geom_rect(xmin=0, xmax=3.5, aes(ymin=-Inf, ymax=Inf), fill='blue', alpha=0.006) +
         geom_rect(xmin=3.5, xmax=5.5, aes(ymin=-Inf, ymax=Inf), fill='grey', alpha=0.006) +
         geom_rect(xmin=5.5, xmax=Inf, aes(ymin=-Inf, ymax=Inf), fill='orange', alpha=0.006) +
         geom_bar(stat='identity', position=position_dodge(), col='black') +
         geom_text(aes(label=deTSSs), position=position_dodge(width=.9), hjust=-.2) +
         labs(title='Differentially expressed CAGE TCs', x='Coefficients', y='Nb. CAGE TCs', caption='Log(FC) >= 1, FDR <= 0.05') +
         coord_flip() + theme_bw() + theme(aspect.ratio=1) +
         scale_y_continuous(expand=c(0.01, 0, 0.15, 0))

# 4f. save
saveRDS(dets, 'data/RDS files/DE_TCs.rds')
saveRDS(TCs, 'data/RDS files/SE_TCs.rds')

