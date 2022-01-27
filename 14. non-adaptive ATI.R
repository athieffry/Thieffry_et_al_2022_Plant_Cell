#### CAGE PTI 2.0: non-adaptive ATI hypothesis
#### Axel Thieffry - August 2021
library(tidyverse)
library(tidylog)
library(magrittr)
library(RColorBrewer)
library(CAGEfightR)
library(vegan)
library(reshape2)
setwd('~/Documents/CAGE 2.0/scripts')
'rename' <- dplyr::rename
'h' <- head
'l' <- length
'at' <- as_tibble


# 1. READ DATA ####
# -----------------
# a. all CAGE TCs
TCs <- readRDS('~/Documents/CAGE 2.0/data/RDS files/SE_TCs.rds')

# b. gene-level CAGE
geneLevel <- readRDS('~/Documents/CAGE 2.0/data/RDS files/SE_genelevel.rds')

# c. ratio of TC expression compared to total gene expression
TC_meanTPM <- assay(TCs, 'TPM') %>% as.data.frame() %>% select(matches('wt_0')) %>% rowMeans() %>% enframe(name='TC_id', value='TC_meanTPM')
gene_meanTPM <- assay(geneLevel, 'TPM') %>% as.data.frame() %>% select(matches('wt_0')) %>% rowMeans() %>% enframe(name='geneID', value='gene_meanTPM')
TC_gene_correspondances <- rowRanges(TCs) %>% as.data.frame() %>% select(geneID) %>% na.omit() %>% rownames_to_column('TC_id') %>% at()

TC_gene_correspondances %>%
    left_join(TC_meanTPM, by='TC_id') %>%
    left_join(gene_meanTPM, by='geneID') %>%
    mutate('ratio'=TC_meanTPM/gene_meanTPM) %>%
    ggplot(aes(x=ratio)) +
           geom_histogram() +
           geom_vline(xintercept=.1)


# 2. COMPUTE TC EXPRESSION RATIO ####
# -----------------------------------
# a. get TC expression
TC_TPM <- assay(TCs, 'TPM') %>%
          as.data.frame() %>%
          rownames_to_column('TC_id') %>%
          melt(id.vars='TC_id', value.name='TPM', variable.name='lib') %>%
          separate('lib', c('genotype', 'time', 'rep'), sep='_') %>%
          group_by(TC_id, genotype, time) %>%
          summarize('mean_TC_TPM'=mean(TPM)) %>%
          unite('sample', c(genotype, time), sep='_') %>%
          pivot_wider(names_from='sample', values_from='mean_TC_TPM') %>%
          ungroup() %>%
          at()

# b. get gene expression
gene_TPM <- assay(geneLevel, 'TPM') %>%
            as.data.frame() %>%
            rownames_to_column('geneID') %>%
            melt(id.vars='geneID', value.name='TPM', variable.name='lib') %>%
            separate('lib', c('genotype', 'time', 'rep'), sep='_') %>%
            group_by(geneID, genotype, time) %>%
            summarize('mean_gene_TPM'=mean(TPM)) %>%
            unite('sample', c(genotype, time), sep='_') %>%
            pivot_wider(names_from='sample', values_from='mean_gene_TPM') %>%
            ungroup() %>%
            at()

# c. TC to gene correspondence
TC_gene_correspondence <- rowRanges(TCs) %>%
                          as.data.frame() %>%
                          select(geneID) %>%
                          filter(!is.na(geneID)) %>%
                          rownames_to_column('TC_id') %>%
                          at()

# d. merge all
TC_TPM %<>% left_join(TC_gene_correspondence, by='TC_id') %>%
            filter(!is.na(geneID)) %>%
            select(TC_id, geneID, everything())

gene_TPM %<>% left_join(TC_gene_correspondence, by='geneID') %>%
              select(geneID, TC_id, everything())

all_merged <- left_join(TC_TPM %>% melt(id.vars=c('TC_id', 'geneID'), variable.name='sample', value.name='TC_TPM'),
                        gene_TPM %>% melt(id.vars=c('TC_id', 'geneID'), variable.name='sample', value.name='gene_TPM'),
                        by=c('TC_id', 'geneID', 'sample')) %>%
              separate('sample', c('genotype', 'time'), sep='_') %>%
              at()


# e. compute Shannon diversity index
diversity <- all_merged %>%
             group_by(geneID, genotype, time) %>%
             summarize('shannon'=vegan::diversity(index='shannon', x=TC_TPM, base=10),
                       'simpson'=vegan::diversity(index='simpson', x=TC_TPM, base=10)) %>%
             ungroup()

# f. add gene expression
diversity %<>% left_join(gene_TPM %>%
                             select(-TC_id) %>%
                             unique() %>%
                             melt(id.vars='geneID', variable.name='lib', value.name='gene_TPM') %>%
                             separate('lib', c('genotype', 'time'), sep='_'),
                         by=c('geneID', 'genotype', 'time'))

# g. genes with multiple TCs
multiTC_genes <- dplyr::count(TC_gene_correspondence, geneID) %>%
                 filter(n!=1) %>%
                 pull(geneID)

# h. plot Simpson index
diversity %>%
    filter(genotype=='wt', time==0) %>%
    ggplot(aes(x=gene_TPM, y=simpson)) +
           geom_point(alpha=.2) +
           scale_x_log10() +
           cowplot::panel_border(color='black') +
           theme(aspect.ratio=1, legend.position='none', axis.line=element_blank(), panel.grid=element_blank()) +
           labs(subtitle='TC diversity of a gene generally decreases with the gene expression level\n(as Xu et al. 2019, Fig 1A)',
                x='Gene expression level (average CAGE TPM, log-scale)',
                y='Simpson index of TC diversity',
                caption='Sample: wild type t=0')

# i. plot Shannon index
diversity %>%
    filter(genotype=='wt', time==0) %>%
    ggplot(aes(x=gene_TPM, y=shannon)) +
           geom_point(alpha=.2) +
           scale_x_log10() +
           cowplot::panel_border(color='black') +
           theme(aspect.ratio=1, legend.position='none', axis.line=element_blank(), panel.grid=element_blank()) +
           labs(subtitle='TC diversity of a gene generally decreases with the gene expression level\n(as Xu et al. 2019, Fig 1C)',
                x='Gene expression level (average CAGE TPM, log-scale)',
                y='Shannon index of TC diversity',
                caption='Sample: wild type t=0')



# 3. FRACTIONAL USAGE OF MAJOR/MINOR TCs ####
# -------------------------------------------
# a. get TSS rank within gene
ranked <- all_merged %>%
          mutate('TC_ratio'=TC_TPM/gene_TPM) %>%
          group_by(geneID, genotype, time) %>%
          mutate('TC_rank'=rank(-TC_ratio, ties.method='random')) %>%
          ungroup()

# b. plot
ranked %>%
    filter(genotype=='wt', time==0, TC_rank < 7) %>%
    mutate('TC_rank'=paste0('TC rank: ', TC_rank)) %>%
    ggplot(aes(x=gene_TPM, y=TC_ratio, col=TC_rank)) +
           geom_point(size=.1, alpha=.2) +
           scale_x_log10() +
           facet_wrap(~TC_rank, ncol=2) +
           theme(aspect.ratio=1, axis.line=element_blank(), panel.grid=element_blank(), legend.position='none') +
           cowplot::panel_border(color='black') +
           scale_color_brewer(palette='Dark2', name='TSS Rank') +
           labs(x='Gene expression level\n(average CAGE TPM, log-scale)',
                y='Fractional usage of TSSs',
                subtitle='Increased fractional usage of the most frequently\nused TSS of a gene and decreased fractional\nusage of each other TSS when gene expression\nlevel rises (as Xu et al. 2019, Fig 2A)',
                caption='Sample: wild type t=0')


