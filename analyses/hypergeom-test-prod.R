# Will Connell
# 2018.10.11

# Depends -----------------------------------------------------------------
devtools::load_all()
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

# calculate
final <- hypergeom_test(cell_meta)
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
plot_pval_distrib(final, p_value_enrich)
plot_pval_distrib(final, p_value_deplete)
