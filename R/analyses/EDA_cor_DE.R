# Will Connell
# 2018.09.16

total_guides <- 282
dummy_frame <- data.frame(guide_cov = c(0:(total_guides-1)))

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
clust10 <- avg_guides %>% 
  filter(louvain == 10) %>% 
  right_join(dummy_frame, by = "guide_cov") %>% 
  select(-louvain) %>% 
  mutate_all(funs(impute.mean(.)))

clust11 <- avg_guides %>% 
  filter(louvain == 11) %>% 
  right_join(dummy_frame, by = "guide_cov") %>% 
  select(-louvain) %>% 
  mutate_all(funs(impute.mean(.)))

covmat <- cor(clust10, clust11)

heatmap(covmat, keep.dendro = FALSE)
