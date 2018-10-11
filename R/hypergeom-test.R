# Will Connell
# 2018.10.10

library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)

# all cell-guide - cluster
path_cell_meta <- "~/ye/scCD4/scanpy/KO_cells.csv"
cell_meta <- read_csv(path_cell_meta) %>% sample_frac(0.5)

# all guide-gene combos
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.txt"
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene) %>% sample_frac(0.5)

# all cluster-gene combos
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.txt"
gene_data_full <- read_tsv(gene_loc_full) %>% 
  select(cluster, gene) %>% sample_frac(0.5)

# intersect DFs
# m, total unique gene-guide comobs
m <- cell_meta %>% 
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>% 
  inner_join(guide_data) %>%
  count(guide, gene) %>% # collapse into guide-gene combos
  mutate(key = paste0(guide, gene))
# n, (total gene-guide combos) - (unique gene-guide-combos)
n <- m %>% 
  mutate(n.hyper = nrow(cell_meta) - n)
# K, size of sample
K <- cell_meta %>% 
  count(louvain)
# q, guide-gene combo w/i cluster
q <- cell_meta %>% 
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>% 
  inner_join(guide_data) %>% 
  count(louvain, guide, gene)

# principle for 
# cluster0-ARID5A.96550280.CCCCGCCGTACCTCTCGTAG_guide-ABCB1
# P(Observed #x cluster-guide-gene or more) 
#  = 1-P(Observed less than #x)
# test 
p.value <- c()
for(cluster in sort(unique(cell_meta$louvain)) ){
  q.clust <- q %>% 
    filter(louvain == cluster) %>% 
    mutate(key = paste0(guide, gene))
  k.clust <- K %>% 
    filter(louvain == cluster)
  m.clust <- m %>% 
    filter(key %in% q.clust$key)
  n.clust <- n %>% 
    filter(key %in% q.clust$key)
  p.value.cluster <- phyper(q = q.clust$n,
                    m = m.clust$n,
                    n = n.clust$n.hyper,
                    k = k.clust$n,
                    lower.tail = FALSE)
  p.value <<- append(p.value, p.value.cluster)
}

final <- as_tibble(cbind(q, p.value))
final.filt <- final %>% 
  mutate(p.adjust = p.adjust(p.value, method = 'fdr')) %>% 
  filter(p.adjust <= 0.01)

# which genes are not shared between clusters?
uq.ge <- list()
for(i in unique(final.filt$louvain)) {
    not.clust.ge <- final.filt %>% 
      filter(louvain != i) %>% 
      select(gene) %>% 
      distinct()
    clust.uq.ge <- final.filt %>% 
      filter(louvain == i) %>% 
      select(gene) %>% 
      anti_join(not.clust.ge)
    colnames(clust.uq.ge) <- sprintf('cluster_%i', i)
    uq.ge <<- append(uq.ge, clust.uq.ge)
}

# which guides are not shared between clusters?
uq.cl <- list()
for(i in unique(final.filt$louvain)) {
  not.clust.guide <- final.filt %>% 
    filter(louvain != i) %>% 
    select(guide) %>% 
    distinct()
  clust.uq.guide <- final.filt %>% 
    filter(louvain == i) %>% 
    select(guide) %>% 
    anti_join(not.clust.guide)
  colnames(clust.uq.guide) <- sprintf('cluster_%i', i)
  uq.cl <<- append(uq.cl, clust.uq.guide)
}
# compared to total # of guides each cluster:
no.unique.guides <- lapply(uq.cl, function(cluster){
  length(cluster)
})
final.filt %>% 
  count(louvain) %>% 
  mutate(diff = nn - unlist(no.unique.guides))

# NEXT
# visualize most important combos in each cluster
# write out csv of guide-gene associated p.vals
# run entire population on cluster (prod scripts)








##### retry rgate method ######
# m, total guide cells
m <- cell_meta %>% 
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>% 
  inner_join(guide_data) %>%
  count(guide) %>% 
  rename(m = n)
# n, total gene cells
n <- cell_meta %>% 
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>% 
  inner_join(guide_data) %>%
  count(gene)
# K, size of sample
K <- cell_meta %>% 
  count(louvain)
# q, guide-gene combo w/i cluster
q <- cell_meta %>% 
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>% 
  inner_join(guide_data) %>% 
  count(guide, gene) %>% 
  rename(q = n)

master <- q %>% 
  inner_join(m) %>% 
  inner_join(n)

# test 
p.value <- c()
for(cluster in sort(unique(cell_meta$louvain)) ){
  k.clust <- K %>% 
    filter(louvain == cluster)
  p.value.cluster <- phyper(q = master$q,
                            m = master$m,
                            n = master$n,
                            k = k.clust$n,
                            lower.tail = FALSE)
  p.value <<- append(p.value, p.value.cluster)
}
# FAIL
final <- cbind(q, p.value)

# correct for multiple comparisons

