# Will Connell
# 2018.09.20

# Monocle Analysis - REPREX


# Depends -----------------------------------------------------------------
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("monocle")
library(monocle)
library(dplyr)


# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../out/%s_monocle", start.date)
dir.create(out.dir)

# Import ------------------------------------------------------------------
data <- readRDS('../../../data_subset_clust_10_11.rds') %>% 
  mutate_if(is.array, as.vector)

expr_matrix <- t(as.matrix(data[11:ncol(data)]))
colnames(expr_matrix) <- data$index
sample_sheet <- data[2:10]
rownames(sample_sheet) <- data$index
gene_annotation <- data.frame(gene_short_name = colnames(data[11:ncol(data)]))
rownames(gene_annotation) <- colnames(data[11:ncol(data)])

# create CellDataSet object
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
# choose a data distribution
# using guassianff() here since data is already normally dist,
# on raw want to use negbinomial.size()
HSMM <- newCellDataSet(as(expr_matrix, 'sparseMatrix'),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = inv.gaussianff()) # must check this inv.gaussian()

# construct single cell trajectories
# use low level of expresison detection - these genes
# account for DE genes already so want to keep them all
# skip cell filtering in general...already preprocessed
HSMM <- detectGenes(HSMM, min_expr = 0.0001)
print('Number cells expressing retained genes summary:')
print(summary(fData(HSMM)$num_cells_expressed))


# METHOD: Construct single cell trajectories ------------------------------
# superset of genes expressed in at least 5% of all cells
fData(HSMM)$use_for_ordering <-
  fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)
# percentage of retained genes
print(sum(fData(HSMM)$use_for_ordering) / length(fData(HSMM)$use_for_ordering))

# visualize importance of PCs by how much
# variation they explain
png(sprintf('%s/pc_variance.png', out.dir), width = 6, height = 6, units = 'in', res = 200)
plot_pc_variance_explained(HSMM, return_all = F)
dev.off()

# reduce the top PCs further using tSNE
HSMM <- reduceDimension(HSMM,
                            max_components = 2,
                            norm_method = 'none',
                            pseudo_expr = 0,
                            num_dim = 4,
                            reduction_method = 'tSNE',
                            verbose = T)

# run density peak clustering to identify the 
# clusters on the 2-D t-SNE space
HSMM <- clusterCells(HSMM, verbose = F)
pData(HSMM)$louvain <- as.factor(pData(HSMM)$louvain)
# facet plot all phenotypically encoded data types
louv_plot <- plot_cell_clusters(HSMM, color_by = 'louvain')
clust_plot <- plot_cell_clusters(HSMM, color_by = 'Cluster')
png(sprintf('%s/subclusters.png', out.dir), width = 8, height = 8, units = 'in', res = 200)
gridExtra::grid.arrange(louv_plot, clust_plot, ncol = 2)
dev.off()

# visualize cell local density (P) vs.
# nearest comparative cell distance (delta)
# for cluster understanding
png(sprintf('%s/rho_delta_thresholds.png', out.dir), width = 6, height = 6, units = 'in', res = 200)
plot_rho_delta(HSMM)
dev.off()

# HSMM <- clusterCells(HSMM,
#                      rho_threshold = 23,
#                      delta_threshold = 23,
#                      skip_rho_sigma = T,
#                      verbose = F)
# louv_plot <- plot_cell_clusters(HSMM, color_by = 'louvain')
# clust_plot <- plot_cell_clusters(HSMM, color_by = 'Cluster')
# gridExtra::grid.arrange(louv_plot, clust_plot, ncol = 2)

# perform DE gene test to extract distinguishing genes
clustering_DEG_genes <- differentialGeneTest(HSMM,
                          fullModelFormulaStr = '~louvain',
                          cores = 1)
# select top significant genes
HSMM_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:150]

HSMM <- setOrderingFilter(HSMM,
                    ordering_genes = HSMM_ordering_genes)

HSMM <- reduceDimension(HSMM, method = 'DDRTree')

HSMM <- orderCells(HSMM)

#HSMM <- orderCells(HSMM, root_state = NULL)

louv_plot <- plot_cell_trajectory(HSMM, color_by = "louvain")
clust_plot <- plot_cell_trajectory(HSMM, color_by = "Cluster")
pseudo_plot <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
png(sprintf('%s/trajectory.png', out.dir), width = 12, height = 8, units = 'in', res = 200)
gridExtra::grid.arrange(louv_plot, clust_plot, pseudo_plot, ncol = 3)
dev.off()

# save object for futher use
saveRDS(HSMM, sprintf("%s/analyzed_obj.RDS", out.dir))


  
