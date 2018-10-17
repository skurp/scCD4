# Will Connell
# 2018.10.11

# Depends -----------------------------------------------------------------
library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)
library(reshape2)
set.seed(55)

# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../../out/%s_hypergeom-test-twoTail", start.date)
dir.create(out.dir)
dir.create("../../out/logs")


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
# K, size of sample, total unique gene-guide combos
K <- cell_meta %>%
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>%
  inner_join(guide_data) %>%
  count(guide, gene) %>% # collapse into guide-gene combos
  mutate(key = paste0(guide, gene))
# m, total number of cells in cluster
m <- cell_meta %>%
  count(louvain)
# n, (total cells NOT in cluster)
n <- m %>%
  mutate(n.hyper = nrow(cell_meta) - n)
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
hypergeom.test <- function(meta) {
  p.value.enrich <- c()
  p.value.deplete <- c()
  p.adj <- c()
  for(cluster in sort(unique(cell_meta$louvain)) ){
    q.clust <- q %>%
      filter(louvain == cluster) %>%
      mutate(key = paste0(guide, gene))
    m.clust <- m %>%
      filter(louvain == cluster)
    k.clust <- K %>%
      filter(key %in% q.clust$key)
    n.clust <- n %>%
      filter(louvain == cluster)
    # test for over representation (enrichment)
    p.value.cluster.enrich <- phyper(q = q.clust$n - 1,
                                     m = m.clust$n,
                                     n = n.clust$n.hyper,
                                     k = k.clust$n,
                                     lower.tail = FALSE)
    # test for under representation (depletion)
    p.value.cluster.deplete <- phyper(q = q.clust$n,
                                      m = m.clust$n,
                                      n = n.clust$n.hyper,
                                      k = k.clust$n,
                                      lower.tail = TRUE)
    # calculating fdr based on ranks in each cluster...
    # validate this thinking
    #p.adj.cluster <- p.adjust(p.value.cluster, method = 'fdr')
    p.value.enrich <- append(p.value.enrich, p.value.cluster.enrich)
    p.value.deplete <- append(p.value.deplete, p.value.cluster.deplete)
    #p.adj <- append(p.adj, p.adj.cluster)
  }
  calc <- data.frame(p.value.enrich = p.value.enrich,
                     p.value.deplete = p.value.deplete)
                     #,p.adjust = p.adj
  as_tibble(cbind(q, calc))
}

# calculate
final <- hypergeom.test(cell_meta)
#final$p.adjust_all <- p.adjust(final$p.value, method = 'fdr')
# write out csv of guide-gene associated p.vals
write_csv(final, sprintf('%s/KO_sigpos_p-vals.csv', out.dir) )

# Q-Q plot of p-values
# png(sprintf("%s/qq-test.png", out.dir), width = 8, height = 9, units = 'in', res = 200)
# plot( x = -log10(ppoints(length(final$p.value))),
#       y = -log10(sort(final$p.value)),
#       xlab= "Expected (-log10)",
#       ylab="Observed (-log10)" )
# abline(0,1,lty=45)
# dev.off()

# histogram of nominal and adjusted p-values
ggplot(final) +
  geom_histogram(aes(x = p.value.enrich), bins = 50) +
  ggtitle("Nominal Enriched p-values") +
  ggsave(sprintf("%s/p-val-enrich.png", out.dir), width = 10, height = 10, units = 'in')
ggplot(final) +
  geom_histogram(aes(x = p.value.deplete), bins = 50) +
  ggtitle("Nominal Depleted p-values") +
  ggsave(sprintf("%s/p-val-deplete.png", out.dir), width = 10, height = 10, units = 'in')
# ggplot(final) +
#   geom_histogram(aes(x = p.adjust), bins = 100) +
#   ggtitle("Adjusted p-values") +
#   ggsave(sprintf("%s/p-adj.png", out.dir), width = 10, height = 10, units = 'in')
# ggplot(final) +
#   geom_histogram(aes(x = p.adjust_all), bins = 100) +
#   ggtitle("Adjusted_all p-values") +
#   ggsave(sprintf("%s/p-adj_all.png", out.dir), width = 10, height = 10, units = 'in')
#
# # filter to examine distribution closely
# filt.final <- final %>% filter(p.adjust < 1)
# ggplot(filt.final) +
#   geom_histogram(aes(x = p.adjust), bins = 100) +
#   xlim(0, 0.5) +
#   ggtitle("Adjusted p-values") +
#   ggsave(sprintf("%s/p-adj-filt.png", out.dir), width = 10, height = 10, units = 'in')
# filt.final <- final %>% filter(p.adjust_all < 1)
# ggplot(filt.final) +
#   geom_histogram(aes(x = p.adjust_all), bins = 100) +
#   xlim(0, 0.5) +
#   ggtitle("Adjusted_all p-values") +
#   ggsave(sprintf("%s/p-adj_all-filt.png", out.dir), width = 10, height = 10, units = 'in')


