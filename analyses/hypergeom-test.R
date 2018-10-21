# Will Connell
# 2018.10.11

# Depends -----------------------------------------------------------------
devtools::load_all()
library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)
library(reshape2)
set.seed(55)

# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../out/%s_hypergeom-test-rprx", start.date)
dir.create(out.dir)
dir.create("../out/logs")


# Input data --------------------------------------------------------------
# all cell-guide - cluster
path_cell_meta <- "~/ye/projects/scanpy/KO_cells.csv"
cell_meta <- read_csv(path_cell_meta) %>% sample_frac(0.05)

# all guide-gene combos
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.txt"
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene) %>% sample_frac(0.05)

# all cluster-gene combos
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.txt"
gene_data_full <- read_tsv(gene_loc_full) %>%
  select(cluster, gene) %>% sample_frac(0.05)

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

# calculate
final <- hypergeom_test(cell_meta, q, m, K, n)
# write out csv of guide-gene associated p.vals
write_csv(final, sprintf('%s/KO_p-vals.csv', out.dir) )

# Q-Q plot of p-values
# png(sprintf("%s/qq-test.png", out.dir), width = 8, height = 9, units = 'in', res = 200)
# plot( x = -log10(ppoints(length(final$p.value))),
#       y = -log10(sort(final$p.value)),
#       xlab= "Expected (-log10)",
#       ylab="Observed (-log10)" )
# abline(0,1,lty=45)
# dev.off()

# histogram of nominal and adjusted p-values
# histogram of nominal and adjusted p-values
melt.final <- final %>%
  melt(measure.vars = c("louvain"),
       id.vars = c("p_value_enrich", "p_value_deplete", "p_adj_enrich", "p_adj_deplete"))
# enriched test
ggplot(melt.final) +
  geom_histogram(aes(x =  p_value_enrich), bins = 100) +
  facet_wrap(vars(value), nrow = 9, scales = "free") +
  ggsave(sprintf("%s/p-val-enrich.png", out.dir), width = 15, height = 30, units = 'in')
# depleted test
ggplot(melt.final) +
  geom_histogram(aes(x =  p_value_deplete), bins = 100) +
  facet_wrap(vars(value), nrow = 9, scales = "free") +
  ggsave(sprintf("%s/p-val-deplete.png", out.dir), width = 15, height = 30, units = 'in')
# FDR enriched test
ggplot(melt.final) +
  geom_histogram(aes(x =  p_adj_enrich), bins = 100) +
  facet_wrap(vars(value), nrow = 9, scales = "free") +
  ggsave(sprintf("%s/q-val-enrich.png", out.dir), width = 15, height = 30, units = 'in')
# FDR depleted test
ggplot(melt.final) +
  geom_histogram(aes(x =  p_adj_deplete), bins = 100) +
  facet_wrap(vars(value), nrow = 9, scales = "free") +
  ggsave(sprintf("%s/q-val-deplete.png", out.dir), width = 15, height = 30, units = 'in')



# Power analysis
# https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html
# determine appropriate effect size (conventional estimations...)


# which genes are not shared between clusters?
uq.ge <- not_shared_feats(filt.final, "gene")
# write out
# write_csv(uq.ge, sprintf('%s/KO_sigpos_uq-gene.csv', out.dir))

# which guides are not shared between clusters?
uq.gu <- not_shared_feats(filt.final, "guide")
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
