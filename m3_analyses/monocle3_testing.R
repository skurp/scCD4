# Will Connell
# 2018.11.08

# Monocle 3 Analysis

# Depends -----------------------------------------------------------------
library(monocle)
library(readr)
library(reticulate)
use_condaenv("monocle3")

# Import ------------------------------------------------------------------
cds <- readRDS("../data/sampled_monocle_obj.RDS")
cds <- updateCDS(cds)
pData(cds)$louvain <- as.factor(pData(cds)$louvain)
pData(cds)$wt <- as.factor(pData(cds)$guide_cov == "0")


# Process Data ------------------------------------------
# Pass TRUE if you want to see progress output on some of Monocle 3's operations
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# project data onto top PCs
cds <- preprocessCDS(cds, num_dim = 20, norm_method = 'none', pseudo_expr = 0,
                     relative_expr = FALSE, scaling = FALSE, verbose = TRUE)
# reduce dimensionality
cds <- reduceDimension(cds, reduction_method = 'UMAP')
# partition cells
cds <- partitionCells(cds)
# learn principle graph
cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
# visualize trajectory
plot_cell_trajectory(cds,
                     color_by = "louvain")

#a helper function to identify the root principal points:
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)

  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
                                           (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
node_ids <- get_correct_root_state(cds, cell_phenotype = "louvain", root_type = "0")
cds <- orderCells(cds, root_pr_nodes = node_ids)
plot_cell_trajectory(cds, color_by = "Pseudotime")

# identfy genes that vary in expression over a trajectory
pr_graph_test <- principalGraphTest(cds, k=3, cores=1)
dplyr::add_rownames(pr_graph_test) %>%
  dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-qval)) %>% head(3)
