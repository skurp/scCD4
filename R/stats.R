# Will Connell
# 2018.10.17


# Statistics --------------------------------------------------------------

# hypergeometric test calculations
hypergeom_test <- function(meta, q, m, K, n) {
  p.value.enrich <- c()
  p.value.deplete <- c()
  p.adj.enrich <- c()
  p.adj.deplete <- c()
  for(cluster in sort(unique(meta[['louvain']])) ){
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
    # calculating fdr based on ranks in each cluster
    p.adj.cluster.enrich <- p.adjust(p.value.cluster.enrich, method = 'fdr')
    p.adj.cluster.deplete <- p.adjust(p.value.cluster.deplete, method = 'fdr')
    # append
    p.value.enrich <- append(p.value.enrich, p.value.cluster.enrich)
    p.value.deplete <- append(p.value.deplete, p.value.cluster.deplete)
    p.adj.enrich <- append(p.adj.enrich, p.adj.cluster.enrich)
    p.adj.deplete <- append(p.adj.deplete, p.adj.cluster.deplete)
  }
  calc <- data.frame(p_value_enrich = p.value.enrich,
                     p_value_deplete = p.value.deplete,
                     p_adj_enrich = p.adj.enrich,
                     p_adj_deplete = p.adj.deplete)
  as_tibble(cbind(q, calc))
}


# Extract features unique to each cluster
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


# Visualizations ----------------------------------------------------------

