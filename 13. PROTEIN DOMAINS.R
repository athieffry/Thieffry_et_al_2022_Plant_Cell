#### ROTEIN DOMAIN ANALYSIS
#### Axel Thieffry - December 2020
set.seed(42)
library(tidyverse)
library(tidylog)
library(ggforce)
library(magrittr)
library(stringr)
library(reshape2)
library(CAGEfightR)
library(RColorBrewer)
library(GenomicRanges)
library(SummarizedExperiment)
library(gprofiler2)
library(eulerr)
library(BiocParallel)
library(patchwork)
library(WriteXLS)
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



# 1. READ INPUTS ####
# -------------------
universe <- readRDS('data/RDS files/universe_geneID.rds')$geneID
myseqinfo <- readRDS('data/RDS files/myseqinfo.rds')
idmapping <- readRDS('data/RDS files/id_mapping.rds')
idmapping2 <- readRDS('data/RDS files/id_mapping2.rds')
intragenicTCs <- readRDS('data/RDS files/SE_intragenicTSSs.rds')
deTCs <- readRDS('data/RDS files/DE_TCs.rds')



# 2. MAP PROTEIN DOMAINS TO TAIR10 ####
# -------------------------------------
# 2a. read TAIR10 protein domains & parse
# source: https://www.arabidopsis.org/download_files/Proteins/Domains/all.domains.txt
prot_domains <- read.csv2('data/ProteinDomains/TAIR10_all.domains', h=F, '\t', strip.white=T, blank.lines.skip=T, stringsAsFactors=F, na.strings=c('', 'NUL', 'NULL')) %>% as_tibble() %>% select(-V2)
    # rename columns
    colnames(prot_domains) <- c('tx_id', 'len', 'app', 'dom_id', 'dom_desc', 'start_dom', 'end_dom', 'evalue', 'date', 'interpro_id', 'interpro_name')
    # clean a bit
    prot_domains %<>% select(-len, -evalue, -date)
    # get geneID from tx_id and re-order
    prot_domains %<>% mutate('geneID'=str_remove(tx_id, '\\.[:digit:]+$')) %>% select(tx_id, geneID, everything())
    # put 'no description' to NA
    prot_domains %<>% mutate('dom_desc'=ifelse(dom_desc=='no description', NA, dom_desc))

# 2b. extract CDS regions from TAIR10 annotation
txdb <- TxDb.Athaliana.BioMart.plantsmart28
cds_regions <- cdsBy(txdb, 'tx', use.names=T)
    # fix seqinfo of CDS regions
    seqlevels(cds_regions) <- seqlevels(myseqinfo)
    seqinfo(cds_regions) <- myseqinfo

# 2c. make IRanges from protein domains
prot_domains_gr <- with(prot_domains, IRanges(start=start_dom * 3, end=end_dom * 3))
mcols(prot_domains_gr) <- prot_domains
names(prot_domains_gr) <- mcols(prot_domains_gr)$dom_id

# 2d. map to genome
on_genome <- pmapFromTranscripts(x=prot_domains_gr, transcripts=cds_regions[mcols(prot_domains_gr)$tx_id])
mcols(on_genome) <- mcols(prot_domains_gr)

# 2e. sort and export
for_export <- on_genome
names(for_export) <- mcols(for_export)$dom_id
for_export <- for_export[width(for_export) > 0]
for_export <- for_export[elementNROWS(for_export) > 0]
export(for_export, 'data/beds/protein_domains_on_genome.bed')



# 3. ANNOTATE TCs WITH DOMAINS ####
# ---------------------------------
# 3a. flatten domains and remove those with negative width
flat_domains <- unlist(on_genome)
flat_domains$tx_id <- names(flat_domains)
flat_domains$geneID <- str_remove(flat_domains$tx_id, '\\.[:digit:]+$')
flat_domains <- subset(flat_domains, width > 0)

# 3b. TCs *in* domain
disruptive_TCs <- rowRanges(intragenicTCs)
hits <- findOverlaps(swapRanges(disruptive_TCs), flat_domains)
disruptive_TCs$in_domain <- extractList(flat_domains$geneID, as(hits, 'List'))
disruptive_TCs$in_domain <- mapply(is.element, disruptive_TCs$geneID, disruptive_TCs$in_domain)

# 3c. TC *downstream* domain
disruptive_TCs$downstream_domain <- follow(swapRanges(disruptive_TCs), flat_domains) # find ANY domain downstream of TC
disruptive_TCs$downstream_domain <- flat_domains$geneID[disruptive_TCs$downstream_domain] # substitute with GENEID of that domain
disruptive_TCs$downstream_domain <- with(disruptive_TCs, ifelse(geneID == downstream_domain, T, F)) # test if protein domain geneID is same as CAGE TC geneID
disruptive_TCs$downstream_domain <- ifelse(is.na(disruptive_TCs$downstream_domain), F, disruptive_TCs$downstream_domain) # change NAs to FALSE

# 3d. protein domain disruption: stats for text
disruptive_TCs %<>% subset(in_domain | downstream_domain)

l(disruptive_TCs) # 454 intragenic TCs
n_distinct(disruptive_TCs$geneID) # in 428 genes
sum(disruptive_TCs$in_domain) # 220 TCs in a domain
subset(disruptive_TCs, in_domain)$geneID %>% n_distinct() # for 216 unique genes
sum(disruptive_TCs$downstream_domain) # 339 TCs downstream a domain
subset(disruptive_TCs, downstream_domain)$geneID %>% n_distinct() # for 317 unique genes

# 3e. make dataframe version
disruptive_TCs_df <- as.data.frame(disruptive_TCs) %>%
                     rownames_to_column('TC_id') %>%
                     as_tibble()

# 3f. proportion of intragenic TCs that are disruptive, by txType category
left_join(disruptive_TCs_df$txType %>% table() %>% enframe(name='txType', value='disruptive_TCs'),
          rowRanges(intragenicTCs)$txType %>% table() %>% enframe(name='txType',value='total_TCs'),
          by='txType') %>%
  subset(disruptive_TCs != 0) %>%
  mutate('percent'=round(disruptive_TCs/total_TCs*100, 2))



# 4. DOMAIN-DISRUPTIVE TCs DE in TIMECOURSE ###
# ---------------------------------------------
# 4a. get DE TCs in timecourse & make wider version of direction
deTCs_timecourse <- deTCs %>% filter(coef %in% c('t10', 't3010', 't30'))

deTCs_timecourse_wide <- deTCs_timecourse %>%
                         select(coef, TC_id, direction) %>%
                         pivot_wider(names_from='coef', values_from='direction')
# 4b. add info to disruptive TCs
disruptive_TCs_df %<>% mutate('TC_DE_in_timecourse'=ifelse(TC_id %in% unique(deTCs_timecourse$TC_id), T, F))
disruptive_TCs_df %<>% left_join(deTCs_timecourse_wide, by='TC_id')

# 4c. stats
# number of disruptive TCs that are DE in timecourse
filter(disruptive_TCs_df, TC_DE_in_timecourse) %>% nrow() # 127 TCs
filter(disruptive_TCs_df, TC_DE_in_timecourse)$geneID %>% n_distinct() # in 125 genes

# 4d. GSEA on genes with induced, domain-disrupting CAGE TC
go_all <- filter(disruptive_TCs_df, TC_DE_in_timecourse) %>%
          pull(geneID) %>%
          unique() %>%
          gost(organism='athaliana', significant=T, user_threshold=0.05, correction_method='fdr', custom_bg=universe, ordered_query=F)

go_all$result %>%
  mutate('p_value_transformed'=-log(p_value, base=10)) %>%
  slice_max(order_by=p_value_transformed, n=10, with_ties=F) %>%
  ggplot(aes(x=10:1, y=p_value_transformed)) +
         geom_col(fill='grey30', alpha=.5) +
         geom_text(aes(label=term_name, y=.1), hjust=0, col='black') +
         coord_flip() +
         cowplot::theme_cowplot() + theme(axis.text.y=element_blank(), aspect.ratio=.75) +
         scale_y_continuous(expand=c(0, 0)) +
         labs(x='Top GO term', y='-log(adj.P.Value, base=10)', title='GSEA',
              subtitle='Genes with domain-disrupting TC induced\nin flg22 treatment (N=125)')

# 4e. export and save
left_join(disruptive_TCs_df, select(idmapping2, -symbol), by='geneID') %>%
  WriteXLS(ExcelFileName='data/excel_outputs/protein_domain_disruptive_TCs.xlsx', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T, verbose=T)

