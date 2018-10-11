# Will Connell
# 2018.10.03

library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)

# cell - guide - cluster
path_cell_meta <- "~/ye/scCD4/scanpy/KO_cells.csv"
cell_meta <- read_csv(path_cell_meta)

# cluster - guide
path_guide_pos <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.sigpos.txt"
path_guide_data <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.vs0.igtb.guide.de.seurate.txt.meta.guide.meta.txt"

guide_pos <- read_tsv(path_guide_pos) %>% 
  rename(guide = cluster) %>% 
  select(guide, gene)
guide_data <- read_tsv(path_guide_data) %>%
  rename(guide = cluster) %>%
  select(guide, gene)

# get DE genes in clusters 3 and 4
gene_loc <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.sigpos.txt"
gene_loc_full <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.txt"
gene_data <- read_tsv(gene_loc) %>% 
  select(cluster, gene)
gene_data_full <- read_tsv(gene_loc_full)

# get sig pos guides present in cluster4
#clust4_guides <- unique(guide_pos$guide)[(unique(guide_pos$guide) %in% unique(clust4$guide_cov))]

# cells split by cluster downsampled to sig pos KO
clusters <- cell_meta %>%
  split(.$louvain) %>% 
  lapply(., function(clust){
    clust <- clust[clust$guide_cov %in% unique(guide_pos$guide), ]
    clust <- clust %>% 
      count(guide_cov) %>% 
      arrange(desc(n))
    clust %>% 
      .[1:50,] %>%  
      mutate(guide_cov = fct_inorder(guide_cov))
  })
# plot
for(i in 1:9){
  clust.num <- i-1
  ggplot(clusters[[i]]) +
    geom_bar(aes(y = n, x = guide_cov), stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(sprintf("Sig-Pos KO Guide Freq Cluster %i", clust.num)) +
    labs(x = "Guide", y = "Freq") +
    coord_flip() +
    ggsave(sprintf('KO_sigpos_freq-%i.png', clust.num), width = 10, height = 20, units = 'in')
}

# take top 50 KO representations from cluster 4,
# extract KO in all other clusters that are NOT found in 4
clust4_ko <- clusters[[5]]
clusters_ko <- clusters
clusters_ko[[5]] <- NULL
clusters_ko <- lapply(clusters_ko, function(clust) {
  clust[!(clust4_ko$guide_cov %in% clust$guide_cov), ]
})

for(i in 1:length(clusters_ko)){
  ggplot(clusters_ko[[i]]) +
    geom_bar(aes(y = n, x = guide_cov), stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(sprintf("Sig-Pos KO Guide Freq NOT in 4; Cluster %s", names(clusters_ko)[i])) +
    labs(x = "Guide", y = "Freq") +
    coord_flip() +
    ggsave(sprintf('KO_sigpos_freq_NOT4-%s.png', names(clusters_ko)[i]), width = 10, height = 20, units = 'in')
}

# get list of all of these unique genes
clusters_ko_unique <- bind_rows(clusters_ko) %>% 
  group_by(guide_cov) %>% 
  summarize(total = sum(n)) %>% 
  arrange(desc(total))

cluster4 <- cell_meta %>%
  filter(louvain == 4) %>% 
  filter(guide_cov %in% clusters_ko_unique$guide_cov[1:15]) %>% 
  count(guide_cov) %>% 
  arrange(desc(n))


cell_meta %>% 
  count(guide_cov) %>% 
  arrange(desc(n))

ctrl <- cell_meta %>% 
  filter(guide_cov == 0) %>% 
  count(louvain)
total <- cell_meta %>% 
  count(louvain)

ctrl$n / total$n


# important
supra <- cell_meta %>% 
  select(index, guide_cov, louvain) %>%
  rename(guide = guide_cov) %>% 
  inner_join(guide_pos)
# fraction cells remaining
length(unique(supra$index)) / length(unique(cell_meta$index))
# compared to fraction of cells that are not ctrls
foo <- cell_meta  %>% filter(guide==0)
1 - (length(unique(foo$index)) / length(unique(cell_meta$index)))
# ...so these are probably mostly ctrls lost when joining by guide
# add genes significantly associated with each cluster
# and subset to genes both sig assoc in cluster and in KO guide
supra <- supra %>% 
  inner_join(gene_data) %>% 
  filter(louvain == cluster)

cluster_genes <- supra %>%
  split(.$louvain) %>%
  lapply(., function(cluster){
    cluster %>% 
      count(gene) %>% 
      arrange(desc(n)) %>% 
      mutate(gene = fct_inorder(gene)) %>% 
      .[1:20,]
  })

for(i in 1:length(cluster_genes)){
  ggplot(cluster_genes[[i]]) +
    geom_bar(aes(y = n, x = gene), stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(sprintf("Sig-Pos Guide & Gene Freq, Cluster %s", names(cluster_genes)[i])) +
    labs(x = "Sig pos gene associated sig pos guide", y = "Freq") +
    coord_flip() +
    ggsave(sprintf('../out/2018-10-09_KO_sigpos_freq/KO_sigpos_freq_guide&gene-%s.png', names(cluster_genes)[i]), width = 10, height = 20, units = 'in')
}


  