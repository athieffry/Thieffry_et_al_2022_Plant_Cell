#### uORF ANALYSIS
#### Axel Thieffry - December 2020
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(reshape2)
library(CAGEfightR)
library(GenomicRanges)
library(BiocParallel)
library(gprofiler2)
library(rtracklayer)
library(readxl)
library(WriteXLS)
library(TxDb.Athaliana.BioMart.plantsmart28)
register(MulticoreParam(workers=6))
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
'count' <- dplyr::count
'l' <- length
'h' <- head


# 0. INPUT DATA ####
# ------------------
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
universe <- readRDS('data/RDS files/universe_geneID.rds')$geneID
deTCs <- readRDS('data/RDS files/DE_TCs.rds')
idmapping2 <- readRDS('data/RDS files/id_mapping2.rds')
intragenic_TCs <- readRDS('data/RDS files/SE_intragenicTSSs.rds')
    # drop txType empty levels
    rowRanges(intragenic_TCs)$txType %<>% droplevels()
txdb <- TxDb.Athaliana.BioMart.plantsmart28
    # fix txdb seqinfo
    seqlevelsStyle(txdb) <- seqlevelsStyle(myseqinfo)
    seqlevels(txdb) <- seqlevels(myseqinfo)
    genome(txdb) <- genome(myseqinfo)
    seqlengths(txdb) <- seqlengths(myseqinfo)



# 1. READ, CLEAN AND PARSE uORF DATA ####
# ---------------------------------------
if(FALSE){
      # 1a. read the set of Arabidopsis thaliana uORFs (from Kurihara et al. 2018)
      uorf <- bind_rows(read_xlsx('data/uORFs_PNAS_2018_Kurihara.xlsx', sheet='WT_Dark'),
                        read_xlsx('data/uORFs_PNAS_2018_Kurihara.xlsx', sheet='WT_Blue'),
                        read_xlsx('data/uORFs_PNAS_2018_Kurihara.xlsx', sheet='hy5_Dark'),
                        read_xlsx('data/uORFs_PNAS_2018_Kurihara.xlsx', sheet='hy5_Blue'))
      # 1b. clean a bit
      uorf %<>% select(-ID:-strand, -annotation, -uTSS_position:-dTSS_TPM)
      uorf %<>% rename('geneID'='gene', 'pep_len'='uORF_peptide_length', 'abs_pos'='absolute_position_uORFs', 'peptide'='uORF_peptide')
      uorf %<>% unique()

      # 1c. de-duplicate uORF gene-wise: keep the gene with most uORFs
      uorf %<>% group_by(geneID) %>%
                filter(n_uORFs==max(n_uORFs)) %>% # keep geneID with the most uORFs
                ungroup() %>%
                filter(!duplicated(geneID)) %>% # in case of ties, keep the first
                arrange(geneID)

      # 1d. verify consistency in provided info & remove erroneous uORF entries
      get_n_peptide  <- function(x) { getElement(x, 'peptide') %>% str_split(', ', simplify=T) %>% length()}  # number of sequences
      get_n_pep_len  <- function(x) { getElement(x, 'pep_len') %>% str_split(', ', simplify=T) %>% length()}  # number of lengths provided
      get_n_pos      <- function(x) { getElement(x, 'abs_pos') %>% str_split(', ', simplify=T) %>% length()}  # number of positions
      get_prov_len   <- function(x) { getElement(x, 'n_uORFs')}                                               # number provided

      counts_df <- tibble('geneID'    = uorf$geneID,
                          'n_peptide' = apply(uorf, 1, get_n_peptide),
                          'n_pep_len' = apply(uorf, 1, get_n_pep_len),
                          'n_pos'     = apply(uorf, 1, get_n_pos),
                          'prov_len'  = apply(uorf, 1, get_prov_len) %>% as.numeric())

      error <- uorf[apply(select(counts_df, -geneID), 1, function(x) var(x) != 0), ] %>% pull(geneID) # still a couple of inconsistent numbers
      uorf %<>% filter(geneID %!in% unique(error)) ; rm(error) # remove inconsistent entries

      # 1e. split uORFs & re-class columns
      uorf %<>% mutate('peptide'=strsplit(peptide, ', '),
                       'pep_len'=strsplit(pep_len, ', '),
                       'abs_pos'=strsplit(abs_pos, ', ')) %>%
                unnest(c(peptide, pep_len, abs_pos))%>%
                mutate('pep_len'=parse_number(pep_len),
                       'abs_pos'=parse_number(abs_pos))

      # 1f. make GenomicRanges object
          # a. calculate uORF length in nucleotides
          uorf %<>% mutate('nt_length'=pep_len * 3)
          # b. add strand and chromosome from geneID
          all_genes_strand <- genes(txdb) %>% as.data.frame() %>% select('geneID'='gene_id', strand, 'chr'=seqnames) %>% as_tibble()
          all(uorf$geneID %in% all_genes_strand$geneID)
          uorf %<>% left_join(all_genes_strand, by='geneID')
          # c. add start & end absolute position (depending on strandness)
          uorf %<>% mutate('start'=ifelse(strand=='+', abs_pos, abs_pos-nt_length+1),
                           'end'=ifelse(strand=='+', start+nt_length-1, abs_pos)) %>%
                    select(-abs_pos)
          # d. do all peptides start with M? YES
          substring(uorf$peptide, 1, 1) %>% table(useNA='ifany')
          # c. make GR
          uorf_gr <- makeGRangesFromDataFrame(uorf, seqnames.field='chr',
                                              start.field='start', end.field='end', strand.field='strand',
                                              keep.extra.columns=T, seqinfo=myseqinfo)
          # d. sort ranges
          uorf_gr %<>% sort()
          # e. export as bed for IGV & genome browser screenshots
          export.bed(unstrand(uorf_gr), 'data/beds/uORFs_7238.bed')
          # f. save RDS
          saveRDS(uorf_gr, file='data/RDS files/uORFs_GR_7238.rds')
          }

uorf_gr <- readRDS('data/RDS files/uORFs_GR_7238.rds')



# 2. POTENTIAL FOR uORF EXCLUSION ####
# ------------------------------------
# 2a. number of genes with an uORF
n_distinct(uorf_gr$geneID) # 3,433 genes

# 2b. length distribution of uORFs
uorf_gr %>%
  as.data.frame() %>%
  ggplot(aes(x=width)) +
         geom_density(fill='grey80') +
         scale_x_log10(expand=c(0, 0)) +
         cowplot::theme_cowplot() + theme(aspect.ratio=.35) +
         scale_y_continuous(expand=c(0, 0)) +
         labs(x='uORF length (bp)', title='uORF length distribution')

# 2c. TCs *in* uORF
uorfloss_TCs <- rowRanges(intragenic_TCs)
hits <- findOverlaps(swapRanges(uorfloss_TCs), uorf_gr)
uorfloss_TCs$in_uORF <- extractList(uorf_gr$geneID, as(hits, 'List'))
uorfloss_TCs$in_uORF <- mapply(is.element, uorfloss_TCs$geneID, uorfloss_TCs$in_uORF)

# 2d. TC *downstream* of uORF
uorfloss_TCs$downstream_uORF <- follow(swapRanges(uorfloss_TCs), uorf_gr) # find nb. ANY TCs downstream of uORF
uorfloss_TCs$downstream_uORF <- uorf_gr$geneID[uorfloss_TCs$downstream_uORF] # add the GENEID of that uORF
uorfloss_TCs$downstream_uORF <- with(uorfloss_TCs, ifelse(geneID == downstream_uORF, T, F)) # conditional: if UPSTREAM GENEID == CURRENT GENEID
uorfloss_TCs$downstream_uORF <- ifelse(is.na(uorfloss_TCs$downstream_uORF), F, uorfloss_TCs$downstream_uORF) # change NAs to FALSE

# 2e. select TCs in or downstream uORF
uorfloss_TCs %<>% subset(in_uORF | downstream_uORF)
l(uorfloss_TCs) # 2,741 TCs are uORF-disruptive (within or downstream uORF)
n_distinct(uorfloss_TCs$geneID) # 2,475 unique genes

# 2f. how many of those are DE in timecourse?
deTCs_timecourse <- deTCs %>% filter(coef %in% c('t10', 't3010', 't30'))
uorfloss_TCs$DE_timecourse <- ifelse(names(uorfloss_TCs) %in% deTCs_timecourse$TC_id, T, F)

subset(uorfloss_TCs, DE_timecourse & (in_uORF | downstream_uORF)) %>% l() # 510 TCs that are DE in the timecourse have potential for uORF loss
subset(uorfloss_TCs, DE_timecourse & (in_uORF | downstream_uORF))$geneID %>% n_distinct() # found in 480 unique genes

# 2g. GO enrichment analysis
genes_with_disrupted_uORF <- subset(uorfloss_TCs, DE_timecourse & (in_uORF | downstream_uORF))$geneID %>% unique()

go_BP <- gost(query=genes_with_disrupted_uORF, organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:BP')
go_MF <- gost(query=genes_with_disrupted_uORF, organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:MF')
go_CC <- gost(query=genes_with_disrupted_uORF, organism='athaliana', ordered_query=F, significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, sources='GO:CC')

gores_BP <- go_BP$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)
gores_MF <- go_MF$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)
gores_CC <- go_CC$result %>% as_tibble() %>% select(p_value, term_name) %>% arrange(p_value) %>% mutate('p_value'=-log(p_value, base=10)) %>% slice_max(order_by=p_value, n=10, with_ties=F)

rbind(mutate(gores_BP, 'set'='BP'),
      mutate(gores_MF, 'set'='MF'),
      mutate(gores_CC, 'set'='CC')) %>%
      ggplot(aes(x=22:1, y=p_value)) +
             geom_col(fill=c(rep('darkgreen', 10), rep('darkred', 2), rep('yellow', 10)), alpha=.5, lwd=.3, col='black') +
             geom_text(aes(label=term_name, y=c(rep(c(.2, .03), each=10), rep(.03, 2))), hjust=0, col='black') +
             facet_wrap(~set, scales='free', ncol=1) +
             coord_flip() + cowplot::theme_cowplot() +
             theme(axis.text.y=element_blank()) +
             scale_y_continuous(expand=c(0, 0)) +
             labs(x='Top10 GO terms', y='-log(adj.P.value, 10)',
                  title='GSEA of genes disrupted uORF', subtitle='Genes with TC in or downstream uORF & DE in timecourse (N=480)')



# 3. OUTPUT RESULTS ####
# ----------------------
# 3a. add DE direction in fgl22 treatement coefficients
deTCs_timecourse_directions <- deTCs_timecourse %>%
                               select(coef, TC_id, direction) %>%
                               pivot_wider(names_from='coef', values_from='direction')

uorfloss_TCs %>%
  as.data.frame() %>%
  rownames_to_column('TC_id') %>%
  filter(DE_timecourse & (in_uORF | downstream_uORF)) %>%
  left_join(deTCs_timecourse_directions, by='TC_id') %>%
  left_join(select(idmapping2, -symbol), by='geneID') %>%
  WriteXLS(ExcelFileName='data/excel_outputs/uORF_disruptive_TCs_DE_timecourse.xlsx', verbose=T, row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)
