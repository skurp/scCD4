# Will Connell
# 2018.10.30

# Guide UpSet charts

# Depends -----------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
devtools::load_all()

guides <- read_tsv("/netapp/home/rgate/tfko_140/nsnp20.raw.sng.guide_sng.meta.guidecluster.enrichment.zsig.txt")

tidy_guides <- apply(guides, 1, function(x) strsplit(x, split = ","))

tidy_feats <- function(tidy_guides, type) {
  df_out <- data.frame(louvain = c(),
                       guide = c())
  for(cluster in tidy_guides){
    feat_type <- cluster[[type]]
    n_feat <- length(feat_type)
    cluster_rep <- as.integer(rep(cluster[['cluster']], n_feat))
    df_tmp <- data.frame(louvain = cluster_rep,
                         guide = feat_type)
    df_out <- rbind(df_out, df_tmp)
  }
  df_out
}

enriched_guides <- tidy_feats(tidy_guides, "enriched")
count(enriched_guides, louvain)
depleted_guides <- tidy_feats(tidy_guides, "depleted")
count(depleted_guides, louvain)

# as compared to all other clusters
not_shared_feats <- function(hypergeom.analysis, feat) {
  uq.feat <- list()
  for( i in sort(unique(hypergeom.analysis[['louvain']]))) {
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
  df <- lapply(uq.feat,`length<-`, max(lengths(uq.feat))) %>%
    as.data.frame()
  colnames(df) <- as.character(sort(unique(hypergeom.analysis[['louvain']])))
  df
}

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
