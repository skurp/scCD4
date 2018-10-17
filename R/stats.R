# Will Connell
# 2018.10.17


# Statistics --------------------------------------------------------------

# hypergeometric test calculations
hypergeom_test <- function(meta) {
  p.value.enrich <- c()
  p.value.deplete <- c()
  p.adj.enrich <- c()
  p.adj.deplete <- c()
  for(cluster in sort(unique(meta['louvain'])) ){
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


# Visualizations ----------------------------------------------------------

# p-value distributions
plot_pval_distrib <- function(pValueDF, p_value_type) {
  plot_list <- list()
  for( cluster in sort(unique(pValues['louvain']))) {
    clust <- final.vals %>%
      filter(louvain == cluster)
    plot <- ggplot(clust) +
      geom_histogram(aes(x = p_value_type), bins = 100) +
      ggtitle(sprintf("Cluster %i", cluster))
    plot_list <- append(plot_list, plot)
  }
  gridExtra::grid.arrange(grobs = plot_list, ncol = 3, nrow = 3) ## display plot
  ggsave(sprintf("%s/%s.png", out.dir, p_value_type),
         arrangeGrob(grobs = plot_list, ncol = 3, nrow = 3),
         width = 15, height = 15, units = 'in')
}
