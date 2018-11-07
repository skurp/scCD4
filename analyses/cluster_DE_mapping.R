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
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.sig.txt"
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
pData(HSMM)$guide_cov <- as.factor(pData(HSMM)$guide_cov)
pData(HSMM)$cluster_guide <- as.factor(paste0(pData(HSMM)$louvain, "-", pData(HSMM)$guide_cov))

tic('DE gene test....')
# perform DE gene test to extract distinguishing genes
clustering_DEG_genes <- differentialGeneTest(HSMM,
                                             fullModelFormulaStr = '~cluster_guide',
                                             cores = 8)
toc()
















# import EN/DE guides
guides <- read_tsv("/netapp/home/rgate/tfko_140/nsnp20.raw.sng.guide_sng.meta.guidecluster.enrichment.zsig.txt")
tidy_guides <- apply(guides, 1, function(x) strsplit(x, split = ","))
enriched_guides <- tidy_feats(tidy_guides, "enriched")
# enriched guides for cluster X
clust1_guide_data <- enriched_guides %>%
  filter(louvain == 1)

# get DE genes associated with each EN guide
# sigpos guide-gene combos
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.sigpos.txt"
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene, fdr) %>%
  mutate(guide = str_extract(guide, "^[:alnum:]+\\.[:digit:]+")) %>%
  filter(guide %in% clust1_guide_data$guide)

png(sprintf("%s/test-traj.png", out.dir), width = 8, height = 8, units = 'in', res = 200)
plot_cell_trajectory(HSMM, color_by = "en.gu.clust1", markers = guide_data$gene)
dev.off()
