# Will Connell
# 2018.09.17

sessionInfo()
# Guide correlation analysis  ---------------------------------------------
# script for production on cluster
# install.packages(c('tidyverse', 'ggridges', 'gplots'), 
#                  repos = 'http://cran.us.r-project.org')
# # retrieve hdf5 tool
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5", suppressUpdates = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(rhdf5)
library(ggridges)
library(gplots)


# Import ------------------------------------------------------------------

# view available strucutre
file_loc <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.km_vb1_default.norm.h5ad"
h5ls(file_loc)

# import meta data
obs <- h5read(file_loc, "/obs") # observations, or rows, including metadata
h5closeAll()
vars <- h5read(file_loc, "/var") # variables, or columns
h5closeAll()
#clust_indices <- which(obs$louvain %in% c(10, 11))
#data <- h5read(file_loc, "/X", index = list(NULL, clust_indices)) #normalized matrix
data <- h5read(file_loc, "/X")
h5closeAll()

# quick summary, check number cells in each cluster
obs %>% 
  count(louvain) %>% 
  print()

# reshape data
#obs_indices <- obs[clust_indices,] #
data <- as.data.frame(t(data))
colnames(data) <- c(vars$index)
# use obs_indices here when testing
data <- cbind(obs, data)

# # get average expression value of entire set
# avg_tot <- mean(as.matrix(data[11:ncol(data)]))
# 
# get avg expression of each guide within each cluster
avg_guides <- data %>%
  select(-index, -donor_cov, -well_cov, -numi_cov, -multiplet_cov,
         -nsnp_cov, -percent_mito, -n_counts) %>%
  group_by(louvain, guide_cov) %>%
  summarize_all(funs(mean)) %>% #avg of each gene expressed within each cluster-guide combo groups
  ungroup()

# weighted representations of each cluster-guide combo
weights <- data %>%
  select(-index, -donor_cov, -well_cov, -numi_cov, -multiplet_cov,
         -nsnp_cov, -percent_mito, -n_counts) %>%
  group_by(louvain, guide_cov) %>%
  tally()

# avg of all expressed genes within each cluster-guide combo group
avg_group_exprs <- avg_guides %>%
  select(-louvain, -guide_cov) %>%
  rowMeans()
# weighted
avg_group_exprs <- (avg_group_exprs * weights$n)
avg_group_exprs <- cbind(avg_guides[1:2], avg_group_exprs)

# # avg cluster expression levels for imputation
# avg_cluster_exprs <- avg_guides %>%
#   group_by(louvain) %>%
#   summarise_all(funs(mean)) %>%
#   ungroup() %>%
#   select(-louvain, -guide_cov) %>%
#   rowMeans() %>%   #impute using mean values
#   #apply(., 1, median) %>%   #impute using median values
#   as.data.frame() %>%
#   mutate(cluster = c(0:11)) %>%
#   rename(avg = '.')
# print('Avg. cluster expression for imputation')
# print(str(avg_cluster_exprs))

# # ridge density plot to visualize guide counts
# data %>% 
#   group_by(louvain, guide_cov) %>% 
#   tally() %>% 
#   add_tally(n) %>% 
#   mutate(normalized = n / nn) %>% 
#   ungroup() %>% 
#   mutate(louvain = as.factor(louvain)) %>%
#   ggplot() +
#   geom_density_ridges(aes(x = guide_cov, y = louvain), bandwidth = 1, panel_scaling = FALSE) + 
#   ggsave('cluster_ridge-density_global-scale.png')
# 
# # ridge density plot to visualize guide counts
# data %>% 
#   group_by(louvain, guide_cov) %>% 
#   tally() %>% 
#   add_tally(n) %>% 
#   mutate(normalized = n / nn) %>% 
#   ungroup() %>% 
#   mutate(louvain = as.factor(louvain)) %>%
#   ggplot() +
#   geom_density_ridges(aes(x = guide_cov, y = louvain), bandwidth = .001, panel_scaling = FALSE) + 
#   ggsave('cluster_ridge-density_global-scale_granular.png')


# Mean guide expression correlation heatmap -------------------------------

# impute NAs with avg_tot
filter_clust <- function(x, df, avg_clust) {
  total_guides <- 282
  new.df <-  data.frame(guide_cov = c(0:(total_guides-1)))
  if(!is.null(avg_clust)) {
    for(i in x) {
      cluster_col <- paste0('cluster_', i)
      new.df <- df %>% 
        filter(louvain == i) %>% 
        left_join(new.df, ., by = c('guide_cov') ) %>% 
        replace_na(list(avg_group_exprs = avg_clust[(avg_clust['cluster'] == i), 'avg'])) %>% 
        select(guide_cov,
               !!cluster_col := avg_group_exprs) %>% 
        left_join(new.df, by = c('guide_cov'))
    }
  } else {
    for(i in x) {
      cluster_col <- paste0('cluster_', i)
      new.df <- df %>% 
        filter(louvain == i) %>% 
        left_join(new.df, ., by = c('guide_cov') ) %>% 
        replace_na(list(avg_group_exprs = 0)) %>% 
        select(guide_cov,
               !!cluster_col := avg_group_exprs) %>% 
        left_join(new.df, by = c('guide_cov'))
    }
  }
  new.df
}

# get nUMI counts in each cluster-guide combo
nUMI_counts <- data %>%
  select(-index, -donor_cov, -well_cov, -multiplet_cov,
         -nsnp_cov, -percent_mito, -n_counts) %>% 
  group_by(louvain, guide_cov) %>% 
  summarize_at(vars(numi_cov), funs(sum)) %>% 
  mutate(numi_cov = as.integer(numi_cov)) %>%
  rename(avg_group_exprs = numi_cov) %>% 
  ungroup()

cluster_col <- paste0('cluster_', c(0:11))
cluster_cor_mat <- filter_clust(0:11, nUMI_counts, NULL) %>% 
  select(!!cluster_col)
cluster_cor_mat_expr <- filter_clust(0:11, avg_group_exprs, NULL) %>% 
  select(!!cluster_col)

# Correllation plots ------------------------------------------------------
library(corrplot)
library(Hmisc)
cor_mat <- cor(cluster_cor_mat, cluster_cor_mat_expr)
#p.mat <- cor.mtest(cluster_cor_mat)$p
#dimnames(p.mat) <- dimnames(cor_mat)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png('cor_combo-nUMI-vs-exprs_impute0_heatmap-pretty.png', width = 6, height = 6, units = 'in', res = 200)
corrplot(cor_mat, method = "color", col = col(200),
         type = "upper", order = "original", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         #p.mat = p.mat, 
         sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE,
         # title and formatting
         title = 'Cluster-Guide nUMI Count \n vs Mean Exprs Correlations',
         mar = c(0,0,4,0))
dev.off()

# cor_mat <- cor(cluster_cor_mat, cluster_cor_mat)
# 
# png('cor_heatmap_weighted-exprs_mean-impute.png', width = 6, height = 6, units = 'in', res = 200)
# heatmap.2(cor_mat, dendrogram = 'none', key = FALSE, density.info = 'none', trace = 'none',
#           Rowv = cluster_col, Colv = FALSE,
#           col= colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(25),
#           cexRow=1, cexCol=1, margins=c(12,8),
#           main = 'Mean Weighted Guide Expression \n Correlation Heatmap')
# dev.off()



