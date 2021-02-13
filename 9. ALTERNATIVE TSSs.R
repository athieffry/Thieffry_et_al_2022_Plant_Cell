#### ALTERNATIVE TSSs
#### Axel Thieffry - January 2021
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(GenomicRanges)
library(rtracklayer)
library(ggalluvial)
library(CAGEfightR)
library(SummarizedExperiment)
library(gprofiler2)
library(BiocParallel)
library(patchwork)
library(readxl)
library(TxDb.Athaliana.BioMart.plantsmart28)
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
universe <- readRDS('data/RDS files/universe_geneID.rds')
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
genelevel <- readRDS('data/RDS files/SE_genelevel.rds')
intragenic_TCs <- readRDS('data/RDS files/SE_intragenicTSSs.rds')
    # remove non-canonical chromosomes
    seqlevels(intragenic_TCs, pruning.mode='coarse') <- paste0('Chr', 1:5)
deTCs <- readRDS('data/RDS files/DE_TCs.rds')
DEGs <- readRDS('data/RDS files/DEGs.rds')
TCs <- readRDS('data/RDS files/SE_TCs.rds')
txdb <- TxDb.Athaliana.BioMart.plantsmart28
    # fix seqinfo
    seqlevelsStyle(txdb) <- seqlevelsStyle(myseqinfo)
    seqlevels(txdb) <- seqlevels(myseqinfo)

# more extensive annotation of genes (AGIs, symbols, description) from the list of UniProt reference proteins in Arabidopsis
if(FALSE){
  idmapping2 <- read.csv('data/peptide_signal_scanning/ids.txt', h=F, sep='\t', quote='', stringsAsFactors=F, strip.white=T, blank.lines.skip=T) %>%
                set_colnames(c('geneID', 'symbol', 'description', 'misc')) %>%
                mutate('geneID'=str_remove(geneID, '\\.[0-9]+')) %>%
                select(geneID, symbol, description) %>%
                as_tibble()
  saveRDS(idmapping2, 'data/RDS files/id_mapping2.rds')
  }

idmapping2 <- readRDS('data/RDS files/id_mapping2.rds')



# 2. NUMBER OF TCs PER GENE ####
# ------------------------------
rowRanges(intragenic_TCs) %>%
  as.data.frame() %>%
  count(geneID) %>%
  ggplot(aes(x=factor(n), fill=n==1)) +
         geom_bar(col='black', lwd=.3) +
         geom_text(stat='count', aes(label=..count..), vjust=-1) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1.3) +
         scale_fill_manual(values=c('black', 'white'), name='Single-TC genes') +
         scale_y_continuous(expand=c(0, 0, 0.1, 0)) +
         labs(x='Number of TCs per gene', y=paste0('Number of genes (N=',rowRanges(intragenic_TCs)$geneID %>% n_distinct(),')'))



# 3. SINGLE/MULTI-TC ANNOTATION IN PC/non-PC GENES ####
# -----------------------------------------------------
# 3a. get PC & non-PC genes & add to data
# source: https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gene_lists/TAIR10_gene_type
gene_types <- read.table('data/TAIR10_gene_type.txt') %>%
              set_colnames(c('geneID', 'type')) %>%
              mutate('type'=ifelse(type=='protein_coding', 'PC-gene', 'non-PC gene')) %>%
              as_tibble()

rowRanges(intragenic_TCs)$gene_type <- gene_types[match(rowRanges(intragenic_TCs)$geneID, gene_types$geneID), ]$type
table(rowRanges(intragenic_TCs)$gene_type, useNA='ifany') # 504 non-PC genes, 20490 PC genes

# 3b. make dataframe for this figure
f2b_df <- rowRanges(intragenic_TCs) %>%
          as.data.frame() %>%
          rownames_to_column('TC_id') %>%
          select(TC_id, txType, geneID, gene_type, score) %>%
          as_tibble()

  # add if major TC (based on pooled TPM expression)
  major_TCs <- f2b_df %>%
               group_by(geneID) %>%
               filter(score==max(score)) %>%
               pull(TC_id) %>%
               unique()

  f2b_df %<>% mutate('major_TC'=ifelse(TC_id %in% major_TCs, TRUE, FALSE))

  # re-categorize annotation into "proximal / promoter / gene-body"
  # basically if TF is not already "proximal or promoter", it's "gene-body".
  f2b_df %<>% mutate('simple_txType'=ifelse(txType %in% c('proximal', 'promoter'), as.character(txType), 'gene body'))
  f2b_df %<>% mutate('simple_txType'=factor(simple_txType, levels=c('proximal', 'promoter', 'gene body')))

  # add if single-TC or multi-TC gene
  single_TC_genes <- count(f2b_df, geneID) %>% filter(n==1) %>% pull(geneID) %>% unique()
  f2b_df %<>% mutate('structure'=ifelse(geneID %in% single_TC_genes, 'single-TC', 'multi-TC'))

  # plot PC-genes first, then non-PC genes, then in Adobe Illustrator
  f2b_df %>%
    filter(gene_type=='PC-gene') %>%
    group_by(major_TC, simple_txType, structure) %>%
    summarize('n'=n()) %>%
    ggplot(aes(x=simple_txType, y=n, fill=major_TC)) +
           geom_col() +
           geom_text(aes(label=n), angle=90, position=position_stack(vjust=.5)) +
           facet_wrap(~structure) +
           scale_y_continuous(expand=c(0, 0, 0.02, 0)) +
           cowplot::theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1))

  # non-PC genes
  f2b_df %>%
    filter(gene_type!='PC-gene') %>%
    group_by(major_TC, simple_txType, structure) %>%
    summarize('n'=n()) %>%
    ggplot(aes(x=simple_txType, y=n, fill=major_TC)) +
           geom_col() +
           geom_text(aes(label=n), angle=90, position=position_stack(vjust=.5)) +
           facet_wrap(~structure) +
           scale_y_continuous(expand=c(0, 0, 0.02, 0)) +
           cowplot::theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1))



# 4. PRIMARY / NON-PRIMARY TCs ####
# ---------------------------------
# 4a. make primary TSS regions (most upstream TAIR10 annotated TSS +/- 100 bp)
primary_regions <- genes(txdb) %>%
                   subset(gene_id %in% unique(rowRanges(intragenic_TCs)$geneID)) %>%
                   promoters(upstream=100, downstream=101) %>%
                   trim()
      # fix seqinfo
      seqlevels(primary_regions) <- seqlevels(myseqinfo)
      genome(primary_regions) <- genome(myseqinfo)
      seqinfo(primary_regions)
      # export for IGV sanity check
      export.bed(unstrand(primary_regions), 'data/beds/primary_regions.bed')

# 4b. find overlaps between intragenic TC (whole tag cluster) & primary regions
hits <- findOverlaps(intragenic_TCs, primary_regions)
    # get GR for all hits
    a <- intragenic_TCs[queryHits(hits)] %>% rowRanges()
    b <- primary_regions[subjectHits(hits)]
    # keep only pairs belonging to the same gene
    in_same_gene <- a$geneID == b$gene_id
    a$primary <- in_same_gene
    # extract TCs overlapping primary regions
    TCs_in_primary_region <- subset(a, primary) %>% names() %>% unique()
    # add info to data
    rowRanges(intragenic_TCs)$primary <- ifelse(names(intragenic_TCs) %in% TCs_in_primary_region, T, F)
    # remove intermediate variables
    rm(a, b, TCs_in_primary_region)

# 4c. stats
rowRanges(intragenic_TCs) %>%
  as.data.frame() %$%
  table('intragenic_txType'=txType, primary)

# 4d. export new intragenic_TCs
saveRDS(intragenic_TCs, 'data/RDS files/SE_intragenicTSSs.rds')



# 5. PRIMARY / NON-PRIMARY TCs DE IN TIMECOURSE ####
# --------------------------------------------------
# 5a. add primary/non-primary info to DE TCs
primary_TCs <- rowRanges(intragenic_TCs) %>%
               subset(primary) %>%
               names() %>%
               unique()

deTCs %<>% mutate('primary'=ifelse(TC_id %in% primary_TCs, T, F))

# 5b. add if TC is intragenic
deTCs %<>% mutate('intragenic'=ifelse(TC_id %in% unique(names(intragenic_TCs)), T, F))

# 5c. plot primary/non-primary DE TCs in timecourse
n_5b <- deTCs %>%
        filter(coef %in% c('t10', 't3010', 't30')) %>%
        filter(intragenic) %>%
        pull(TC_id) %>%
        n_distinct()

deTCs %>%
  filter(coef %in% c('t10', 't3010', 't30')) %>%
  filter(intragenic) %>%
  select(coef, primary) %>%
  mutate('primary'=ifelse(primary, 'primary TC', 'alternative TC')) %>%
  group_by(coef, primary) %>%
  summarize('n'=n()) %>%
  mutate('coef'=factor(coef, levels=c('t30', 't3010', 't10'))) %>%
  ggplot(aes(x=coef, y=n, fill=primary)) +
         geom_col(col='black', lwd=.3) +
         geom_text(aes(label=n), position=position_stack(vjust=.5), col=rep(c('black', 'white'), 3)) +
         coord_flip() + cowplot::theme_cowplot() + theme(aspect.ratio=.3) +
         scale_y_continuous(expand=c(0.01, 0, 0, 0.01)) +
         labs(y=paste0('Intragenic DE TCs in timecourse (N=', n_5b, ' unique TCs)'),
              x='Comparison', title='Primary/non-primary TCs DE in timecourse') +
         scale_fill_manual(values=c('white', 'grey50'))

  # proportion of non-primary DE TCs per coef (for text)
  deTCs %>%
    filter(coef %in% c('t10', 't3010', 't30')) %>%
    filter(intragenic) %>%
    group_by(coef) %>%
    summarize('primary_%'=mean(primary)*100,
              'non_primary_%'=mean(!primary)*100)

# 5d. how many genes have at least ONE alternative DE TC in the timecourse?
genes_with_alternative_DE_TC_timecourse <- deTCs %>%
                                           filter(coef %in% c('t10', 't3010', 't30')) %>%
                                           filter(!primary & intragenic) %>%
                                           pull(geneID) %>%
                                           unique()

l(genes_with_alternative_DE_TC_timecourse) # 619

  # proportion of timecourse DEGs that also have an alternative DE TC in timecourse?
  degs_timecourse <- DEGs %>% filter(coef %in% c('t10', 't3010', 't30')) %>% pull(geneID) %>% unique()

  table(degs_timecourse %in% genes_with_alternative_DE_TC_timecourse)
  mean(degs_timecourse %in% genes_with_alternative_DE_TC_timecourse)

# 5e.as above, but split by coefficient
# genes with at least one alternative DE TC in timecourse
deTCs %>%
  filter(coef %in% c('t10', 't3010', 't30')) %>%
  filter(!primary & intragenic) %>%
  select(coef, geneID) %>%
  group_by(coef) %>%
  mutate('coef'=factor(coef, levels=c('t30', 't3010', 't10'))) %>%
  summarize('n_unique_genes_with_alt_DE_TC'=n_distinct(geneID)) %>%
  ggplot(aes(x=coef, y=n_unique_genes_with_alt_DE_TC)) +
         geom_col(fill='grey50', lwd=.3, col='black') +
         geom_text(aes(label=n_unique_genes_with_alt_DE_TC), position=position_stack(vjust=.5), col='white') +
         cowplot::theme_cowplot() + coord_flip() +
         theme(aspect.ratio=.3) + scale_y_continuous(expand=c(0.01, 0, 0, 0.01)) +
         labs(x='Comparision',
              y=paste0('Genes with alternative DE TC in time course\n(N=',l(genes_with_alternative_DE_TC_timecourse),' unique genes)'))

# 5f. number of genes by number of alternative DE TCs in timecourse
deTCs %>%
  filter(coef %in% c('t10', 't3010', 't30')) %>%
  filter(!primary & intragenic) %>%
  select(coef, geneID) %>%
  group_by(coef, geneID) %>%
  summarize('n'=n()) %>%
  filter(n > 1) %>%
  mutate('coef'=factor(coef, levels=c('t30', 't3010', 't10'))) %>%
  ungroup() %>%
  group_by(coef) %>%
  summarize('n_genes_2ormore'=n_distinct(geneID)) %>%
  ggplot(aes(x=coef, y=n_genes_2ormore)) +
         geom_col(fill='grey50', lwd=.3, col='black') +
         geom_text(aes(label=n_genes_2ormore), position=position_stack(vjust=.5), col='white') +
         cowplot::theme_cowplot() + coord_flip() +
         theme(aspect.ratio=.3) + scale_y_continuous(expand=c(0.01, 0, 0, 0.01)) +
         labs(x='Comparision',
              y=paste0('Genes with more than 1 alternative TC\ninduced in time course (N=36 unique genes)'))

# 5g. get genes with 2 or more alternative DE TCs in timecourse
deTCs %>%
  filter(coef %in% c('t10', 't3010', 't30')) %>%
  filter(!primary & intragenic) %>%
  select(coef, geneID) %>%
  group_by(coef, geneID) %>%
  summarize('n'=n()) %>%
  ungroup() %>%
  filter(n > 1) %>%
  select(geneID) %>%
  unique() %>%
  left_join(idmapping2, by='geneID')



# 6. GENES WITH TC SWITCHES (ONE TC UP, ONE TC DOWN) ####
# -------------------------------------------------------
# 6a. get genes with TC switches during timecourse
    t10_switches <- deTCs %>%
                    filter(coef=='t10' & intragenic) %>%
                    select(geneID, TC_id, txType, direction) %>%
                    arrange(geneID)

    t10_tokeep <- t10_switches %>%
                  group_by(geneID) %>%
                  summarize('switch'=n_distinct(direction)==2) %>%
                  filter(switch) %>%
                  pull(geneID)

    t10_switches %<>% filter(geneID %in% t10_tokeep)

    t3010_switches <- deTCs %>%
                    filter(coef=='t3010' & intragenic) %>%
                    select(geneID, TC_id, txType, direction) %>%
                    arrange(geneID)

    t3010_tokeep <- t3010_switches %>%
                  group_by(geneID) %>%
                  summarize('switch'=n_distinct(direction)==2) %>%
                  filter(switch) %>%
                  pull(geneID)

    t3010_switches %<>% filter(geneID %in% t3010_tokeep)

    t30_switches <- deTCs %>%
                    filter(coef=='t30' & intragenic) %>%
                    select(geneID, TC_id, txType, direction) %>%
                    arrange(geneID)

    t30_tokeep <- t30_switches %>%
                  group_by(geneID) %>%
                  summarize('switch'=n_distinct(direction)==2) %>%
                  filter(switch) %>%
                  pull(geneID)

    t30_switches %<>% filter(geneID %in% t30_tokeep)

    # remove intermediates & merge results
    rm(t10_tokeep, t3010_tokeep, t30_tokeep)

    switches <- rbind(t10_switches %>% mutate('coef'='t10'),
                      t3010_switches %>% mutate('coef'='t3010'),
                      t30_switches %>% mutate('coef'='t30'))

    n_distinct(switches$geneID) # 21 genes have a TC switch in timecourse

# 6b. annotation of TC switches (Figure 2G)
switches %>%
  mutate('coef'=factor(coef, levels=c('t10', 't3010', 't30'))) %>%
ggplot(aes(x=txType, fill=direction)) +
       geom_bar(position=position_dodge(), lwd=.3, col='black') +
       facet_wrap(~coef) +
       cowplot::theme_cowplot() +
       scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks()) +
       scale_fill_brewer(palette='Set1', direction=-1) +
       coord_flip()

# 6c. export results
switches %>%
  left_join(deTCs %>%
              filter(coef %in% c('t10', 't3010', 't30')) %>%
              select(coef, TC_id, logFC, adj.P.Val, symbol, primary),
            by=c('TC_id', 'coef')) %>%
  left_join(select(idmapping2, -symbol), by='geneID') %>%
  WriteXLS::WriteXLS('data/excel_outputs/TC_switches_timecourse.xlsx', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)

