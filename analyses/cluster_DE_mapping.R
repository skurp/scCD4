# Will Connell
# 2018.10.19


# Depends -----------------------------------------------------------------
library(tictoc)
library(monocle)
library(dplyr)
library(readr)
set.seed(55)

# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../out/%s_TEST", start.date)
dir.create(out.dir)
dir.create("../out/logs")

# Input data --------------------------------------------------------------
HSMM <- readRDS('../data/sampled_monocle_obj.RDS')
# all cluster-gene combos
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.txt"
gene_data_full <- read_tsv(gene_loc_full)

# METHOD: Construct single cell trajectories ------------------------------
tic('Gene abnd cell filtering.....')

# use low level of expresison detection - these genes
# account for DE genes already so want to keep them all
# skip cell filtering in general...already preprocessed
HSMM <- detectGenes(HSMM, min_expr = 0.0001)
print('Number cells expressing retained genes summary:')
print(summary(fData(HSMM)$num_cells_expressed))

# superset of genes expressed in at least 5% of all cells
fData(HSMM)$use_for_ordering <-
  fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)
# percentage of retained genes
print('Percentage of retained genes:')
print(sum(fData(HSMM)$use_for_ordering) / length(fData(HSMM)$use_for_ordering))
toc()

# visualize
pData(HSMM)$louvain <- as.factor(pData(HSMM)$louvain)
pData(HSMM)$guide <- as.factor(pData(HSMM)$louvain)

cells_to_be_tested <- pData(HSMM) %>%
  filter(louvain == 8)
genes_to_be_tested <- gene_data_full %>%
  filter(cluster == 2)

cds_subset <- HSMM[genes_to_be_tested$gene,]

tic('DE gene test....')
# perform DE gene test to extract distinguishing genes
clustering_DEG_genes <- differentialGeneTest(cds_subset,
                                             fullModelFormulaStr = '~louvain',
                                             cores = 4)
toc()
























# all cell-guide - cluster
path_cell_meta <- "~/ye/projects/scanpy/KO_cells.csv"
cell_meta <- read_csv(path_cell_meta) %>% sample_frac(.1)

# all guide-gene combos
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.txt"
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene, fdr) %>% sample_frac(.1)

# all cluster-gene combos
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.txt"
gene_data_full <- read_tsv(gene_loc_full) %>%
  select(cluster, gene, fdr) %>% sample_frac(.1)

clust8 <- cell_meta %>%
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>%
  inner_join(guide_data) %>%
  filter(louvain == 8) %>%
  filter(guide == "TCF25.89900721.AGCTCTGCTCACACTCAGGG_guide")

foo <- readRDS("../data/data_subset_clust_10_11.rds")
clust10 <- foo %>%
  filter(louvain == 10) %>%
  filter(guide_cov == 0)
all <- foo %>%
  filter(guide_cov == 0)


all <- cell_meta %>%
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>%
  inner_join(guide_data) %>%
  #filter(louvain == 8) %>%
  filter(guide == "TCF25.89900721.AGCTCTGCTCACACTCAGGG_guide")

length(unique(all$gene))



clustA.guideB <- foo %>%
  filter(louvain == 10 & guide_cov == 4)
clustA.WT <- foo %>%
  filter(louvain == 10 & guide_cov == 0)





melt.final <- gene_data_full %>%
  melt(measure.vars = c("p_val1"),
       id.vars = c("cluster"))

# enriched test
ggplot(melt.final) +
  geom_histogram(aes(x =  value), bins = 100) +
  facet_wrap(vars(cluster), nrow = 9, scales = "free") +
  ggsave(sprintf("%s/fdr-distrib.png", out.dir), width = 15, height = 30, units = 'in')


ggplot(gene_data_full) +
  geom_histogram(aes(x =  pvalue), bins = 100) +
  #facet_wrap(vars(cluster), nrow = 9, scales = "free") +
  ggsave(sprintf("%s/pval-TEST.png", out.dir), width = 15, height = 30, units = 'in')
