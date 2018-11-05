# Will Connell
# 2018.10.30

# Guide UpSet charts

# Depends -----------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
devtools::load_all()

guides <- read_tsv("/netapp/home/rgate/tfko_140/nsnp20.raw.sng.guide_sng.meta.guidecluster.enrichment.zsig.txt")

tidy_guides <- apply(guides, 1, function(x) strsplit(x, split = ","))

enriched_guides <- tidy_feats(tidy_guides, "enriched")
count(enriched_guides, louvain)
depleted_guides <- tidy_feats(tidy_guides, "depleted")
count(depleted_guides, louvain)

uq.en <- not_shared_feats(enriched_guides, "guide") %>%
  gather(key = "louvain", value = "guide") %>%
  na.omit() %>%
  mutate(louvain = as.integer(louvain)) %>%
  arrange(louvain)
uq.de <- not_shared_feats(depleted_guides, "guide") %>%
  gather(key = "louvain", value = "guide") %>%
  na.omit() %>%
  mutate(louvain = as.integer(louvain)) %>%
  arrange(louvain)


# Test DEG ordering (trajectory) based on uniquely ENRCICHED guides ------------
# meta file
path_cell_meta <- "~/ye/projects/scanpy/KO_cells.csv"
cell_meta <- read_csv(path_cell_meta) %>%
  rename(guide = guide_cov) %>%
  mutate(guide = str_extract(guide, "^[:alnum:]+\\.[:digit:]+")) %>%
  mutate(guide = if_else(is.na(guide), '0', guide))

# sigpos guide-gene combos
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.sigpos.txt"
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene, fdr) %>%
  mutate(guide = str_extract(guide, "^[:alnum:]+\\.[:digit:]+"))

# uniquely enriched guides for cluster X and DE gene associations
clust1_guide_data <- guide_data %>%
  filter(guide %in% (filter(uq.en, louvain == 1)[['guide']]))
# enriched guides for cluster X
clust1_guide_data <- enriched_guides %>%
  filter(louvain == 1)

# count of cells in clsuter X containing enriched guides
cell_meta %>%
  filter(louvain == 1 & guide %in% clust1_guide_data$guide) %>%
  dim()

