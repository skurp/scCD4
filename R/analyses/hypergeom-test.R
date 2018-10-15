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
out.dir <- sprintf("../../out/%s_hypergeom-test-rprx", start.date)
dir.create(out.dir)
dir.create("../../out/logs")


# Input data --------------------------------------------------------------
# all cell-guide - cluster
path_cell_meta <- "~/ye/projects/scanpy/KO_cells.csv"
cell_meta <- read_csv(path_cell_meta) %>% sample_frac(0.1)

# all guide-gene combos
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.txt"
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene) %>% sample_frac(0.1)

# all cluster-gene combos
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.txt"
gene_data_full <- read_tsv(gene_loc_full) %>%
  select(cluster, gene) %>% sample_frac(0.1)

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
  p.value <- c()
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
    p.value.cluster <- phyper(q = q.clust$n,
                              m = m.clust$n,
                              n = n.clust$n.hyper,
                              k = k.clust$n,
                              lower.tail = FALSE)
    # calculating fdr based on ranks in each cluster...
    # validate this thinking
    p.adj.cluster <- p.adjust(p.value.cluster, method = 'fdr')
    p.value <- append(p.value, p.value.cluster)
    p.adj <- append(p.adj, p.adj.cluster)
  }
  calc <- data.frame(p.value = p.value,
                     p.adjust = p.adj)
  as_tibble(cbind(q, calc))
}

# calculate
final <- hypergeom.test(cell_meta)
final$p.adjust_all <- p.adjust(final$p.value, method = 'fdr')
# write out csv of guide-gene associated p.vals
# write_csv(final, sprintf('%s/KO_sigpos_p-vals.csv', out.dir) )

# Q-Q plot of p-values
png(sprintf("%s/qq-test.png", out.dir), width = 8, height = 9, units = 'in', res = 200)
plot( x = -log10(ppoints(length(final$p.value))),
      y = -log10(sort(final$p.value)),
      xlab= "Expected (-log10)",
      ylab="Observed (-log10)" )
abline(0,1,lty=45)
dev.off()

# Power analysis
# https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html
# determine appropriate effect size (conventional estimations...)
#

# take filter set
final.filt <- final %>%
  filter(p.adjust2 <= 0.01)

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
    colnames(clust.uq.feat) <- sprintf('%i', i)
    uq.feat <- append(uq.feat, clust.uq.feat)
  }
  # corece to dataframe
  lapply(uq.feat,`length<-`, max(lengths(uq.feat))) %>%
    as.data.frame()
}

# which genes are not shared between clusters?
uq.ge <- not_shared_feats(final.filt, "gene")
# write out
# write_csv(uq.ge, sprintf('%s/KO_sigpos_uq-gene.csv', out.dir))

# which guides are not shared between clusters?
uq.cl <- not_shared_feats(final.filt, "guide")
# write out
# write_csv(uq.cl, sprintf('%s/KO_sigpos_uq-guide.csv', out.dir))

print("Number of guides in each cluster:")
final.filt %>%
  count(louvain, guide) %>%
  print()
print("Number of genes in each cluster:")
final.filt %>%
  count(louvain) %>%
  print()

# visualize no. genes associate w/ each guide in each cluster
print(paste("Number of distinct genes associate with all guides:",
            n_distinct(guide_data$gene)) ) # for reference
# why is this one guide missing??
print(
  cell_meta$guide_cov[which(!(unique(cell_meta$guide_cov) %in% unique(guide_data$guide)))]
)

clust_guide <- final.filt %>%
  count(louvain, guide)
ggplot(clust_guide) +
  geom_bar(aes(x = guide, y = nn), stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Sig-Pos KO Guide-Gene Freq") +
  labs(x = "Guide", y = "No. Genes") +
  coord_flip() +
  facet_wrap(vars(louvain), nrow = 1) +
  ggsave(sprintf('%s/KO_sigpos_freq.png', out.dir), width = 30, height = 20, units = 'in')
