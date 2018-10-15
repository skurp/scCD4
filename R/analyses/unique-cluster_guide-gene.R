# Will Connell
# 2018.10.11

# Unique guide-genes between cluster

# Depends -----------------------------------------------------------------
library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)
library(reshape2)
set.seed(55)


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
write_csv(uq.ge, sprintf('%s/KO_sigpos_uq-gene.csv', out.dir))

# which guides are not shared between clusters?
uq.cl <- not_shared_feats(final.filt, "guide")
# write out
write_csv(uq.cl, sprintf('%s/KO_sigpos_uq-guide.csv', out.dir))

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
