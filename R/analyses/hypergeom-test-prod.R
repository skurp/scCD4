# Will Connell
# 2018.10.11

# Depends -----------------------------------------------------------------
library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)

# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../out/%s_hypergeom-test", start.date)
dir.create(out.dir)


# Input data --------------------------------------------------------------
# all cell-guide - cluster
path_cell_meta <- "~/ye/projects/scanpy/KO_cells.csv"
cell_meta <- read_csv(path_cell_meta) %>% sample_frac(1)

# all guide-gene combos
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.txt"
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene) %>% sample_frac(1)

# all cluster-gene combos
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.txt"
gene_data_full <- read_tsv(gene_loc_full) %>%
  select(cluster, gene) %>% sample_frac(1)

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
final <- final %>%
  mutate(p.adjust = p.adjust(p.value, method = 'fdr'))
# write out csv of guide-gene associated p.vals
write_csv(final, sprintf('../../out/%s/KO_sigpos_p-vals.csv', out.dir) )

# take filter set
final.filt <- final %>%
  filter(p.adjust <= 0.01)

# function
# to extract features unique to each cluster
# as compared to all other clusters
not_shared_feats <- function(hypergeom.analysis, feat) {
  uq.feat <- list()
  for( i in sort(unique(hypergeom.analysis$louvain))) {
    not.clust.feat <- hypergeom.analysis %>%
      filter(louvain != i) %>%
      select(feat) %>%
      distinct()
    clust.uq.feat <- hypergeom.analysis %>%
      filter(louvain == i) %>%
      select(feat) %>%
      distinct() %>%
      anti_join(not.clust.feat)
    colnames(clust.uq.feat) <- sprintf('cluster_%i', i)
    uq.feat <- append(uq.feat, clust.uq.feat)
  }
  # corece to dataframe
  lapply(uq.feat,`length<-`, max(lengths(uq.feat))) %>%
    as.data.frame()
}

# which genes are not shared between clusters?
uq.ge <- not_shared_feats(final.filt, "gene")
# write out
write_csv(uq.ge, sprintf('../../out/%s/KO_sigpos_uq-gene.csv', out.dir))

# which guides are not shared between clusters?
uq.cl <- not_shared_feats(final.filt, "guide")
# write out
write_csv(uq.cl, sprintf('../../out/%s/KO_sigpos_uq-guide.csv', out.dir))

# compared to total # of guides each cluster:
no.uq.guides <- uq.cl %>%
  t() %>%
  as_tibble() %>%
  tibble::rownames_to_column("louvain") %>%
  rename(guide = V1) %>%
  count(louvain)
print("Number of guides that are not unique to cluster:")
final.filt %>%
  count(guide) %>%
  mutate(diff = nn - (no.uq.guides$n)) %>%
  cbind(., no.uq.guides$louvain) %>%
  print()

# visualize no. genes associate w/ each guide in each cluster
print(paste("Number of distinct genes associate with all guides:",
            n_distinct(guide_data$gene)) ) # for reference
# why is this one guide missing??
cell_meta$guide_cov[which(!(unique(cell_meta$guide_cov) %in% unique(guide_data$guide)))]

clust_guide <- final.filt %>%
  count(louvain, guide)
ggplot(clust_guide) +
  geom_bar(aes(x = guide, y = nn), stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Sig-Pos KO Guide-Gene Freq") +
  labs(x = "Guide", y = "No. Genes") +
  coord_flip() +
  facet_wrap(vars(louvain), nrow = 1) +
  ggsave(sprintf('../../out/%s/KO_sigpos_freq.png', out.dir), width = 20, height = 10, units = 'in')

