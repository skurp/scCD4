# Will Connell
# 2018.09.20

# Monocle Analysis - REPREX


# Depends -----------------------------------------------------------------
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("monocle")

library(rhdf5)
library(monocle)
library(dplyr)
library(tictoc)
library(readr)


# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../out/%s_monocle-mat2-rprx", start.date)
dir.create(out.dir)

args <- commandArgs(trailingOnly = TRUE)

print('Number of cores:')
print(detectCores())

# Import ------------------------------------------------------------------
tic('Importing & Preparing Data.....')

# view available strucutre
# file_loc <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.km_vb1_default.norm.h5ad"
#  ----- 2018.09.28: UPDATE NEW MATRIX -----------
file_loc <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.h5ad"
h5ls(file_loc)

# import meta data
path_cell_meta <- "~/ye/projects/scanpy/KO_cells.csv"
obs <- read_csv(path_cell_meta) # observations, or rows, including metadata
vars <- h5read(file_loc, "/var") # variables, or columns

# subset a sample to effectively
clust_indices <- sample(1:nrow(obs), as.numeric(args[1]))
h5closeAll()
data <- h5read(file_loc, "/X", index = list(NULL, clust_indices)) #normalized matrix
obs_indices <- obs[clust_indices,]

data <- data %>%
  t() %>%
  as.data.frame()
colnames(data) <- c(vars$index)
data <- cbind(obs_indices, data) %>%
  mutate_if(is.array, as.vector)

print("Dimensions of data:")
print(dim(data))


# Prepare for CellDataSet object ---------------------------------------------
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

toc()


# METHOD: Construct single cell trajectories ------------------------------
tic('Gene abnd cell filtering.....')

# use low level of expresison detection - these genes
# account for DE genes already so want to keep them all
# skip cell filtering in general...already preprocessed
HSMM <- detectGenes(HSMM, min_expr = 0.0001)
print('Number cells expressing retained genes summary:')
print(summary(fData(HSMM)$num_cells_expressed))

# superset of genes expressed in at least 5% of all cells
fData(HSMM)$use_for_ordering <-
  fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)
# percentage of retained genes
print('Percentage of retained genes:')
print(sum(fData(HSMM)$use_for_ordering) / length(fData(HSMM)$use_for_ordering))
toc()

# visualize importance of PCs by how much
# variation they explain
png(sprintf('%s/pc_variance.png', out.dir), width = 6, height = 6, units = 'in', res = 200)
plot_pc_variance_explained(HSMM, return_all = F)
dev.off()

tic('Calc PCs and reduce via tSNE...')
# reduce the top PCs further using tSNE
HSMM <- reduceDimension(HSMM,
                        max_components = 5,
                        norm_method = 'none',
                        pseudo_expr = 0,
                        num_dim = 4,
                        reduction_method = 'tSNE',
                        cores = round(detectCores()*0.9),
                        verbose = T)
toc()

tic('Density peak clustering.....')
# run density peak clustering to identify the
# clusters on the 2-D t-SNE space
HSMM <- clusterCells(HSMM,
                     verbose = T,
                     cores = round(detectCores()*0.9),
                     skip_rho_sigma = T,
                     num_clusters = 9)
                     # delta_threshold = 5,
                     # rho_threshold = 620,
                     # method = "louvain",
                     # k = round(sqrt(dim(HSMM)[[2]])),
                     # louvain_iter = 3)
toc()
# visualize
pData(HSMM)$louvain <- as.factor(pData(HSMM)$louvain)
pData(HSMM)$guide_cov <- as.factor(pData(HSMM)$guide_cov)
pData(HSMM)$wt <- as.factor(pData(HSMM)$guide_cov == "0")
# facet plot all phenotypically encoded data types
louv_plot <- plot_cell_clusters(HSMM, color_by = 'louvain')
clust_plot <- plot_cell_clusters(HSMM, color_by = 'Cluster')
#guide_plot <- plot_cell_clusters(HSMM, color_by = 'guide_cov')
wt_plot <-  plot_cell_clusters(HSMM, color_by = 'wt')
png(sprintf('%s/subclusters.png', out.dir), width = 30, height = 20, units = 'in', res = 200)
gridExtra::grid.arrange(clust_plot, louv_plot, wt_plot, ncol = 3)
dev.off()

# visualize cell local density (P) vs.
# nearest comparative cell distance (delta)
# for cluster understanding
# Can only employ this plot when using
# Density Peak clustering
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

tic('DE gene test....')
# perform DE gene test to extract distinguishing genes
clustering_DEG_genes <- differentialGeneTest(HSMM,
                                             fullModelFormulaStr = '~louvain',
                                             #reducedModelFormulaStr = '~louvain',
                                             cores = round(detectCores()*0.9))
toc()

# select top significant genes
HSMM_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]

HSMM <- setOrderingFilter(HSMM,
                          ordering_genes = HSMM_ordering_genes)

tic('Reduce dimensions via DDRTree.....')
HSMM <- reduceDimension(HSMM, method = 'DDRTree', cores = round(detectCores()*0.9))
toc()

tic('Order cells according to learned trajectory of pseudotime.....')
HSMM <- orderCells(HSMM)
toc()

# visualize
louv_plot <- plot_cell_trajectory(HSMM, color_by = "louvain")
clust_plot <- plot_cell_trajectory(HSMM, color_by = "Cluster")
state_plot <- plot_cell_trajectory(HSMM, color_by = "State")
pseudo_plot <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
#guide_plot <- plot_cell_trajectory(HSMM, color_by = "guide_cov")
wt_plot <- plot_cell_trajectory(HSMM, color_by = "wt")
png(sprintf('%s/trajectory.png', out.dir), width = 30, height = 20, units = 'in', res = 200)
gridExtra::grid.arrange(louv_plot, clust_plot, state_plot,
                        pseudo_plot, wt_plot,
                        ncol = 3, nrow = 2)
dev.off()

# save object for futher use
saveRDS(HSMM, sprintf("%s/analyzed_obj.RDS", out.dir))

