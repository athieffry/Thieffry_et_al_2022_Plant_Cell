#### CAGE Gene-level Differential Expression Analysis
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
design <- readRDS('data/RDS files/design.rds')
contrasts <- readRDS('data/RDS files/contrasts.rds')
genelevel <- readRDS('data/RDS files/SE_genelevel.rds')



# 2. LIMMA DE ANALYSIS ####
# -------------------------
# 2a. run differential expression analysis
dge <- DGEList(counts=assay(genelevel, 'counts'),
               samples=colData,
               lib.size=colData(genelevel)$totalTags) %>%
       calcNormFactors(method='TMM')
v <- voom(dge, design=design, plot=F)
fit <- lmFit(v, design=design)
fit <- contrasts.fit(fit, contrasts)
eb <- eBayes(fit, robust=T)
res <- decideTests(eb, adjust.method='BH', lfc=1, p.value=0.05)

# 2b. get numeric matrix of normalized expression values on the log2 scale
assay(genelevel, 'E') <- v$E

# 2c. extract results for all coefficients
coefs <- colnames(res)
names(coefs) <- coefs

degs <- map(coefs, topTable, fit=eb, number=Inf, adjust.method='BH', p.value=0.05, lfc=1) %>%
  lapply(function(x) rownames_to_column(x, 'geneID')) %>%
  bind_rows(.id='coef') %>%
  mutate('direction'=ifelse(logFC > 0, 'up', 'down')) %>%
  as_tibble()

all_degs <- map(coefs, topTable, fit=eb, number=Inf, adjust.method='BH') %>%
  lapply(function(x) rownames_to_column(x, 'geneID')) %>%
  bind_rows(.id='coef') %>%
  as_tibble()



# 3. DEA POST-PROCESSING ####
# ---------------------------
# 3a. make coefficient factor level
coef_levels <- factor(c('t10', 't3010', 't30', 'hen2', 'rrp4', 'hen2.t10', 'hen2.3010', 'hen2.t30', 'rrp4.t30'))

# 3b. add effect (flg22, exosome, interaction)
degs %<>% mutate('effect'=case_when(.$coef %in% c('hen2', 'rrp4') ~ 'exosome',
                                    .$coef %in% c('t10', 't30', 't3010') ~ 'flg22',
                                    .$coef %in% c('rrp4.t30') ~ 'interaction',
                                    TRUE ~ 'ERROR'))
all_degs %<>% mutate('effect'=case_when(.$coef %in% c('hen2', 'rrp4') ~ 'exosome',
                                        .$coef %in% c('t10', 't30', 't3010') ~ 'flg22',
                                        .$coef %in% c('rrp4.t30', 'hen2.3010', 'hen2.t10', 'hen2.t30') ~ 'interaction',
                                        TRUE ~ 'ERROR'))

    # sanity check
    with(degs, table(coef, effect))
    with(all_degs, table(coef, effect))

# 3c. add DE summary results to TCs SE
metadata(genelevel)$res <- as.matrix(res)

# 3d. plot number of up/down CAGE TCs by coefficient
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
         labs(title='DEGs', x='Coefficients', y='Nb. CAGE TCs', caption='Log(FC) >= 1, FDR <= 0.05') +
         coord_flip() + theme_bw() + theme(aspect.ratio=1) +
         scale_y_continuous(expand=c(0.01, 0, 0.15, 0))

# 3e. save
saveRDS(degs, 'data/RDS files/DEGs.rds')
saveRDS(genelevel, 'data/RDS files/SE_genelevel.rds')
