# Will Connell
# 2018.10.16

# Plotting p-values correctly

final.vals <- read_csv("../../out/2018-10-16_hypergeom-test-twoTail/KO_sigpos_p-vals.csv")

# enriched
for( i in sort(unique(final.vals$louvain))) {
  clust <- final.vals %>%
    filter(louvain == i)
  plot <- ggplot(clust) +
            geom_histogram(aes(x = p.value.enrich), bins = 100)

  png(sprintf("%s/p-val-enrich-clust-%i.png", out.dir, i), width = 10, height = 10, units = 'in', res = 200)
  print(plot)
  dev.off()
}

# depleted
for( i in sort(unique(final.vals$louvain))) {
  clust <- final.vals %>%
    filter(louvain == i)
  plot <- ggplot(clust) +
    geom_histogram(aes(x = p.value.deplete), bins = 100)

  png(sprintf("%s/p-val-deplete-clust-%i.png", out.dir, i), width = 10, height = 10, units = 'in', res = 200)
  print(plot)
  dev.off()
}

final.vals %>%
  #filter(louvain == 0) %>%
  count(louvain, guide)
