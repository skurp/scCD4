library(rhdf5)
library(monocle)
library(dplyr)
library(tictoc)
library(readr)
library(stringr)
devtools::load_all()

# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../out/%s_monocle-mat2-clust1-rprx", start.date)
dir.create(out.dir)

# load analyzed monocle object
HSMM <- readRDS("../out/2018-11-07_monocle-mat2-clust1-rprx/analyzed_obj.RDS")

#all cell-guide - cluster
path_cell_meta <- "~/ye/projects/scanpy/KO_cells.csv"
# NA is WT or 0 guides
cell_meta <- read_csv(path_cell_meta) %>%
  mutate(guide_cov = str_extract(guide_cov, "^[:alnum:]+\\.[:digit:]+"))
top_guides <- cell_meta %>%
  filter(louvain == 1) %>%
  count(guide_cov) %>%
  arrange(desc(n))

# vectorize top frequency guides in cluster 1
pData(HSMM)$top.gu.clust1 <- replace(as.vector(pData(HSMM)$guide_cov),
                                    !(as.vector(pData(HSMM)$guide_cov)%in% top_guides$guide_cov[2:4]),
                                    "non-top-guides")
pData(HSMM)$top.gu.clust1 <- factor(pData(HSMM)$top.gu.clust1,
                                     levels = c("non-top-guides", top_guides$guide_cov[2:4]))

# rearrange "under" points first, "over" points last
pheno <- pData(HSMM)
pheno$rownames <- rownames(pData(HSMM))
pheno <- arrange(pheno, top.gu.clust1)
# pheno$alpha <-  ifelse(pheno$top.gu.clust1 == "non-top-guides", 1, 2)
pData(HSMM) <- pheno
rownames(pData(HSMM)) <- pheno$rownames

# rearrange columns in reducedDim
# to match pData upon merging
cds_reorder <- reducedDimS(HSMM)
cds_reorder <- cds_reorder[,pheno$rownames]
HSMM@reducedDimS <- cds_reorder

# visualize
louv_plot <- plot_cell_trajectory(HSMM, color_by = "louvain")
#clust_plot <- plot_cell_trajectory(HSMM, color_by = "Cluster")
state_plot <- plot_cell_trajectory(HSMM, color_by = "State")
pseudo_plot <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
guide_plot <- plot_cell_trajectory(HSMM, color_by = "top.gu.clust1", size = pheno$alpha)
wt_plot <- plot_cell_trajectory(HSMM, color_by = "wt")
png(sprintf('%s/trajectory.png', out.dir), width = 30, height = 20, units = 'in', res = 200)
gridExtra::grid.arrange(louv_plot, state_plot, #clust_plot,
                        pseudo_plot, wt_plot, guide_plot,
                        ncol = 3, nrow = 2)
dev.off()


# one off plot
test <- pData(HSMM) %>%
  count(donor_cov, guide_cov) %>%
  arrange(desc(n))

foo <- priority_feature_reorder(cds = HSMM, pheno_feature = "donor_cov", priority_vector = test$donor_cov[1])

png(sprintf('%s/trajectory.png', out.dir), width = 8, height = 8, units = 'in', res = 200)
plot_cell_trajectory(foo, color_by = "priority_reorder")
dev.off()



