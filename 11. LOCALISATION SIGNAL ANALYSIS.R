#### LOCALISATION SIGNAL (LS) ANALYSIS
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
library(SummarizedExperiment)
library(BiocParallel)
library(patchwork)
library(limma)
library(edgeR)
library(gprofiler2)
library(WriteXLS)
library(eulerr)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(BSgenome.Athaliana.TAIR.TAIR9)
register(MulticoreParam(workers=6))
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
'count' <- dplyr::count
'l' <- length
'h' <- head


# READ INPUTS ####
# -------------------
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
universe <- readRDS('data/RDS files/universe_geneID.rds')$geneID
idmapping <- readRDS('data/RDS files/id_mapping.rds')
idmapping2 <- readRDS('data/RDS files/id_mapping2.rds')

dtus_timecourse <- readRDS('data/RDS files/DTUs_timecourse.rds')

intragenicTCs <- readRDS('data/RDS files/SE_intragenicTSSs.rds')
    # remove intragenic TCs not on canonical chromosomes
    seqlevels(intragenicTCs, pruning.mode='coarse') <- paste0('Chr', 1:5)

TCs <- readRDS('data/RDS files/SE_TCs.rds')
deTCs <- readRDS('data/RDS files/DE_TCs.rds')

txdb <- TxDb.Athaliana.BioMart.plantsmart28
    # fix seqinfo
    seqlevelsStyle(txdb) <- seqlevelsStyle(myseqinfo)
    seqlevels(txdb) <- seqlevels(myseqinfo)




# 1. DOWNLOAD ALL PROTEIN SEQUENCES FROM ARABIDOPSIS.ORG AND RUN SIGNALP AND TARGETP ####
# ---------------------------------------------------------------------------------------
# run scanning with: SignalP-5.0: http://www.cbs.dtu.dk/services/SignalP and TargetP-2.0: http://www.cbs.dtu.dk/services/TargetP
# nohup signalp -fasta proteins.faa -format 'short' -org 'euk' -plot 'none' -stdout -batch 30000 > signalp.results.txt 2> sig.err > sig.out &
# nohup targetp -fasta proteins.faa -format 'short' -org 'pl' -plot 'none' -stdout -batch 30000 > targetp.results.txt 2> tar.err > tar.out &
# result file: remove 1st line and uncomment column headers before reading into R (just easier)



# 2. COMPILE PREDICTED LOCALISATION SIGNALS ####
# ----------------------------------------------
# 2a. read signalP and targetP results
signalp <- read.table('data/peptide_signal_scanning/signalp.results.txt', h=T, sep='\t') %>%
           rename('tx_id'='ID', 'prediction'='Prediction') %>%
           mutate('prediction'=str_remove(prediction, '\\(Sec/SPI\\)')) %>%
           select(tx_id, prediction, CS.Position) %>%
           as_tibble()

targetp <- read.table('data/peptide_signal_scanning/targetp.results.txt', h=T, sep='\t') %>%
           rename('tx_id'='ID', 'prediction'='Prediction') %>%
           select(tx_id, prediction, CS.Position) %>%
           as_tibble()

# 2b. keep predicted results only
targetp %<>% filter(prediction != 'noTP')
signalp %<>% filter(prediction != 'OTHER')

# 2c. venn diagram for genes with a signal peptide
data.frame('tx_id'=unique(c(signalp$tx_id, targetp$tx_id))) %>%
  mutate('signalP'=ifelse(tx_id %in% signalp$tx_id, 1, 0),
         'targetP'=ifelse(tx_id %in% targetp$tx_id, 1, 0)) %>%
  column_to_rownames('tx_id') %>%
  euler(shape='ellipse') %>% plot(quantities=T, main='All TXs with predicted localisation signal')

# 2d. given venn diagram, add signalP-only results to targetp
locsig <- rbind(targetp, signalp[signalp$tx_id %!in% targetp$tx_id, ])
rm(signalp, targetp)
    # sanity check: none is duplicated
    duplicated(locsig$tx_id) %>% sum()

# 2e. remove entries without cleavage site (N=6) and not on canonical chromosome
locsig %<>% filter(str_detect(CS.Position, '\\?', negate=T))
canonical_chr_idx <- substring(locsig$tx_id, first=3, last=3) %in% 1:5
locsig <- locsig[canonical_chr_idx, ]
rm(canonical_chr_idx)

# 2f. map locsig onto genomic coordinates
    # get CDS coordinates for TAIR10 transcript
    cds <- cdsBy(txdb, by='tx', use.names=T)
    # keep only CDS for TXs with localisation signal
    cds <- cds[names(cds) %in% locsig$tx_id]
    # extract CS position and multiply by 3 for AA to nt conversion
    locsig$CS.Position %<>% str_remove('CS pos: ') %>%
                            str_remove('CS pos luTP: ') %>%
                            str_remove('. Pr: .*') %>%
                            str_remove('-.*') %>%
                            as.integer() %>%
                            multiply_by(3)
    # make locsig IR
    locsig_ir <- with(locsig, IRanges(start=1, end=CS.Position))
    mcols(locsig_ir) <- locsig
    # map locsig to transcript space
    locsig_gr <- pmapFromTranscripts(x=locsig_ir, transcripts=cds[mcols(locsig_ir)$tx_id])
    # remove 0-width regions
    locsig_gr <- locsig_gr[width(locsig_gr) > 0]
    # remove empty elements
    locsig_gr <- locsig_gr[elementNROWS(locsig_gr) > 0]
    # make flat version for use in alternative TSS analysis
    locsig_flat <- locsig_gr %>% range() %>% unlist()
    mcols(locsig_flat) <- data.frame('tx_id'=names(locsig_flat),
                                     'geneID'=str_remove(names(locsig_flat), '\\.[:digit:]+$')) %>%
                          left_join(select(locsig, -CS.Position), by='tx_id')
    # clean and export BED files
    rm(locsig_ir, locsig)
    export.bed(locsig_gr, 'data/beds/localisation_signals.bed')
    export.bed(locsig_flat, 'data/beds/localisation_signals_flat.bed')



# 3. SUMMARY STATS ####
# ---------------------
# 3a. number of predicted localisation signals
l(locsig_flat) # 6985 transcripts
# 3b. distribution of locsig types
table(locsig_flat$prediction) # 1470 cTP, 131 luTP, 1102 mTP, 4282 SP
# 3c. overall GO enrichment
go_all_locsig_BP <- gost(query=locsig_flat$geneID, organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:BP')
go_all_locsig_MF <- gost(query=locsig_flat$geneID, organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:MF')
go_all_locsig_CC <- gost(query=locsig_flat$geneID, organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:CC')
gores_all_locsig_BP <- go_all_locsig_BP$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)
gores_all_locsig_MF <- go_all_locsig_MF$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)
gores_all_locsig_CC <- go_all_locsig_CC$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)

rbind(mutate(gores_all_locsig_BP, 'set'='BP'),
      mutate(gores_all_locsig_MF, 'set'='MF'),
      mutate(gores_all_locsig_CC, 'set'='CC')) %>%
      ggplot(aes(x=30:1, y=p_value)) +
             geom_col(fill=rep(c('darkgreen', 'darkred', 'yellow'), each=10), alpha=.5, lwd=.3, col='black') +
             geom_text(aes(label=term_name, y=rep(c(2, 4, 10), each=10)), hjust=0, col='black') +
             facet_wrap(~set, scales='free') +
             coord_flip() + cowplot::theme_cowplot() +
             theme(axis.text.y=element_blank()) +
             scale_y_continuous(expand=c(0, 0)) +
             labs(x='Top10 GO terms', y='-log(adj.P.value, 10)',
                  title='GSEA of all genes with predicted localisation signal', subtitle='N=6895')



# 4. CAGE TCs WITHIN, DOWNSTREAM, OR UPSTREAM LOC SIG ####
# --------------------------------------------------------
# 4a. CAGE TCs *in* LS
LS_TCs <- rowRanges(intragenicTCs)
hits <- findOverlaps(swapRanges(LS_TCs), locsig_flat)
LS_TCs$in_LS <- extractList(locsig_flat$geneID, as(hits, 'List'))
LS_TCs$in_LS <- mapply(is.element, LS_TCs$geneID, LS_TCs$in_LS)

subset(LS_TCs, in_LS) %>% names() %>% n_distinct() # 184 CAGE TCs *in* LS
subset(LS_TCs, in_LS)$geneID %>% n_distinct() # in 181 genes

# 4b. TC *upstream* of LS
LS_TCs$upstream_LS <- precede(swapRanges(LS_TCs), locsig_flat)
LS_TCs$upstream_LS <- locsig_flat$geneID[LS_TCs$upstream_LS]
LS_TCs$upstream_LS <- with(LS_TCs, ifelse(geneID == upstream_LS, T, F))
LS_TCs$upstream_LS <- ifelse(is.na(LS_TCs$upstream_LS), F, LS_TCs$upstream_LS)

subset(LS_TCs, upstream_LS) %>% names() %>% n_distinct() # 4778 CAGE TCs *downstream* LS
subset(LS_TCs, upstream_LS)$geneID %>% n_distinct() # in 4577 genes

# 4c. TC *downstream* of LS
LS_TCs$downstream_LS <- follow(swapRanges(LS_TCs), locsig_flat) # find nb. any LS downstream of TC
LS_TCs$downstream_LS <- locsig_flat$geneID[LS_TCs$downstream_LS] # add the GENEID of that TC-disruptive LS
LS_TCs$downstream_LS <- with(LS_TCs, ifelse(geneID == downstream_LS, T, F)) # conditional: if UPSTREAM GENEID == CURRENT GENEID
LS_TCs$downstream_LS <- ifelse(is.na(LS_TCs$downstream_LS), F, LS_TCs$downstream_LS) # change NAs to FALSE

subset(LS_TCs, downstream_LS) %>% names() %>% n_distinct() # 109 CAGE TCs *downstream* LS
subset(LS_TCs, downstream_LS)$geneID %>% n_distinct() # in 105 genes

    # in or downstream LS stats for text
    subset(LS_TCs, in_LS | downstream_LS) %>% names() %>% n_distinct() # 293 disrupting CAGE TCs
    subset(LS_TCs, in_LS | downstream_LS)$geneID %>% n_distinct() # in 284 genes

# 4d. GSEA of genes with CAGE TC *in* or *downstream* LS
go_inOrDownstream_LS_BP <- gost(query=unique(subset(LS_TCs, in_LS | downstream_LS)$geneID), organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:BP')
go_inOrDownstream_LS_MF <- gost(query=unique(subset(LS_TCs, in_LS | downstream_LS)$geneID), organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:MF')
go_inOrDownstream_LS_CC <- gost(query=unique(subset(LS_TCs, in_LS | downstream_LS)$geneID), organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:CC')

gores_inOrDownstream_LS_BP <- go_inOrDownstream_LS_BP$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)
gores_inOrDownstream_LS_MF <- go_inOrDownstream_LS_MF$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)
gores_inOrDownstream_LS_CC <- go_inOrDownstream_LS_CC$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)

rbind(mutate(gores_inOrDownstream_LS_BP, 'set'='BP'),
      mutate(gores_inOrDownstream_LS_MF, 'set'='MF'),
      mutate(gores_inOrDownstream_LS_CC, 'set'='CC')) %>%
      ggplot(aes(x=30:1, y=p_value)) +
             geom_col(fill=rep(c('darkgreen', 'darkred', 'yellow'), each=10), alpha=.5, lwd=.3, col='black') +
             geom_text(aes(label=term_name, y=rep(c(0.2, .4, .4), each=10)), hjust=0, col='black') +
             facet_wrap(~set, scales='free') +
             coord_flip() + cowplot::theme_cowplot() +
             theme(axis.text.y=element_blank()) +
             scale_y_continuous(expand=c(0, 0)) +
             labs(x='Top10 GO terms', y='-log(adj.P.value, 10)',
                  title='GSEA of genes TC in or downstream localisation signal', subtitle='N=284')

# 4e. add whether TC is DE in timecourse or not
deTCs_timecourse <- deTCs %>% filter(coef %in% c('t10', 't3010', 't30'))
LS_TCs$TC_DE_timecourse <- names(LS_TCs) %in% unique(deTCs_timecourse$TC_id)

subset(LS_TCs, upstream_LS & TC_DE_timecourse) %>% l() # 922 TCs upstream a LS are DE in timecourse
subset(LS_TCs, upstream_LS & TC_DE_timecourse)$geneID %>% n_distinct() # in 887 genes

subset(LS_TCs, in_LS & TC_DE_timecourse) %>% l() # 44 TCs in a LS are DE in timecourse
subset(LS_TCs, in_LS & TC_DE_timecourse)$geneID %>% n_distinct() # in 43 genes

subset(LS_TCs, downstream_LS & TC_DE_timecourse) %>% l() # 34 TCs downstream a LS are DE in timecourse
subset(LS_TCs, downstream_LS & TC_DE_timecourse)$geneID %>% n_distinct() # in 33 genes

subset(LS_TCs, (in_LS | downstream_LS) & TC_DE_timecourse) %>% l() # 78 TCs downstream a LS are DE in timecourse
subset(LS_TCs, (in_LS | downstream_LS) & TC_DE_timecourse)$geneID %>% n_distinct() # in 76 genes


# 4f. add whether TC is a DTU in timecourse
LS_TCs$DTU_timecourse <- names(LS_TCs) %in% unique(dtus_timecourse$TC_id)

subset(LS_TCs, upstream_LS & DTU_timecourse) %>% l() # 70 TCs are upsteram a LS and timecourse DTU
subset(LS_TCs, upstream_LS & DTU_timecourse)$geneID %>% n_distinct() # in 47 genes

subset(LS_TCs, in_LS & DTU_timecourse) %>% l() # 8 TCs are in a LS and timecourse DTU
subset(LS_TCs, in_LS & DTU_timecourse)$geneID %>% n_distinct() # in 8 genes

subset(LS_TCs, downstream_LS & DTU_timecourse) %>% l() # 18 TCs are downstream a LS and timecourse DTU
subset(LS_TCs, downstream_LS & DTU_timecourse)$geneID %>% n_distinct() # in 18 genes

subset(LS_TCs, (in_LS | downstream_LS) & DTU_timecourse) %>% l() # 26 TCs are downstream a LS and timecourse DTU
subset(LS_TCs, (in_LS | downstream_LS) & DTU_timecourse)$geneID %>% n_distinct() # in 26 genes

# 4g. get genes with timecourse DTU in or downstream LS
subset(LS_TCs, DTU_timecourse & (in_LS | downstream_LS)) %>%
  as.data.frame() %>%
  rownames_to_column('TC_id') %>%
  left_join(idmapping2, by='geneID') %>%
  select(TC_id, txType, primary, upstream_LS, in_LS, downstream_LS, TC_DE_timecourse, DTU_timecourse, geneID, symbol.x, symbol.y, description) %>%
  WriteXLS('data/excel_outputs/LS_DTUs_inOrDownstream.xlsx', row.names=F, col.names=T, BoldHeaderRow=T, AdjWidth=T)



# 5. EXPORT ALL ####
# ------------------
# 5a. make master data frame
LS_TCs_df <- LS_TCs %>%
             as.data.frame() %>%
             filter(in_LS | downstream_LS) %>%
             rownames_to_column('TC_id') %>%
             select(-width, -score, -thick.start:-txID, -txType_extended:-shape_mean, -geneID_anti, -symbol_anti, -composition, -gene_type) %>%
             as_tibble()
    # sanity check
    n_distinct(LS_TCs_df$geneID) # 284

# 5b. add timecourse differential expression results
DE_results_df <- deTCs_timecourse %>%
                 select(coef, TC_id, direction) %>%
                 pivot_wider(names_from='coef', values_from='direction') %>%
                 select(TC_id, 'DE_t10'=t10, 'DE_t3010'=t3010, 'DE_t30'=t30)

LS_TCs_df %<>% left_join(DE_results_df, by='TC_id')

# 5c. add timecourse DTU results
DTU_results_df <- dtus_timecourse %>%
                  select(coef, TC_id, direction) %>%
                  pivot_wider(names_from='coef', values_from='direction') %>%
                  select(TC_id, 'DTU_t10'=t10, 'DTU_t3010'=t3010, 'DTU_t30'=t30)

LS_TCs_df %<>% left_join(DTU_results_df, by='TC_id')

WriteXLS(LS_TCs_df, ExcelFileName='data/excel_outputs/disruptive_LS_results.xlsx', verbose=T, row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)
