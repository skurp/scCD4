# Will Connell
# 2018.09.25

# Correlation matrix b/w cluster-guide
# expression averages

library(dplyr)
library(purrr)
library(corrplot)


# Import ------------------------------------------------------------------
data <- readRDS("../../../data_subset_clust_10_11.rds")
dummy_df <- data.frame(guide_cov = c(0:281))

# get avg expression of each guide within each cluster 
avg_guides <- data %>%
  select(-index, -donor_cov, -well_cov, -numi_cov, -multiplet_cov,
         -nsnp_cov, -percent_mito, -n_counts) %>% 
  group_by(louvain, guide_cov) %>%
  summarize_all(funs(mean)) %>% #avg of each gene expressed within each groups
  ungroup()

# fake data
fake_9 <- data %>%
  select(-index, -donor_cov, -well_cov, -numi_cov, -multiplet_cov,
         -nsnp_cov, -percent_mito, -n_counts) %>% 
  group_by(louvain, guide_cov) %>%
  summarize_all(funs(mean)) %>% #avg of each gene expressed within each groups
  ungroup() %>% 
  filter(louvain == 10) %>% 
  mutate(louvain = louvain - 1) %>% 
  mutate_all(funs(.^2))





filter_clust <- function(x, df) {
  total_guides <- 282
  dummy_df <-  data.frame(guide_cov = c(0:(total_guides-1)))
  for(i in x) {
    if(!exists('new.df')){
      new.df <- df %>% 
        filter(louvain == i) %>% 
        left_join(dummy_df, ., by = c('guide_cov') ) %>% 
        replace_na(list(louvain = i)) %>% 
        replace(is.na(.), 0)
    } else {
      new.df <- df %>% 
        filter(louvain == i) %>% 
        left_join(dummy_df, ., by = c('guide_cov') ) %>% 
        replace_na(list(louvain = i)) %>% 
        replace(is.na(.), 0) %>% 
        bind_rows(new.df)
    }
  }
  new.df
}


foo <- filter_clust(10:11, avg_guides)

bar <- foo %>% 
  select(-guide_cov) %>% 
  group_by(louvain) %>% 
  cor()

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png('TEST.png', width = 20, height = 20, units = 'in', res = 200)
corrplot(bar, method = "color", col = col(200),
         type = "upper", order = "original", number.cex = .7, number.digits = 5,
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         #p.mat = p.mat, 
         sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE,
         # title and formatting
         title = 'Cluster-Guide Mean \n Expression Correlations',
         mar = c(0,0,4,0))
dev.off()
  