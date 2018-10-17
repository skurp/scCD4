# Will Connell
# 2018.09.27

# Test scfind package to compare gene
# lists to HCA database

# get genes pseudtotemporally associated 
# with each cluster


# HACK TEST ---------------------------------------------------------------
library(scfind)
library(SingleCellExperiment)

# get gene expression matrix associated with cluster 11
HSMM <- readRDS("../out/2018-09-24_monocle/analyzed_obj.RDS")

# valid_cells <- row.names(subset(pData(HSMM),
#                                 louvain == 11))
# HSMM <- HSMM[,valid_cells]

diff_test_res <- differentialGeneTest(HSMM,
                                      fullModelFormulaStr = "~sm.ns(louvain)")
#diff_test_res[,c("gene_short_name", "pval", "qval")]
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
str(sig_gene_names)

png("tst.png", width = 8, height = 18, units = 'in', res = 200)
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 2,
                        cores = core_set,
                        show_rownames = T)
dev.off()
gene_list <- c("BCL2A1", "VPS13B", "SLC37A3", "NELL2", "SERPINB9", "FGFBP2",
               "MAML2", "ID3", "PLAG1", "ARSK", "KLHL28")

# SCFIND ------------------------------------------------------------------
data <- readRDS('../../../data_subset_clust_10_11.rds') %>% 
  mutate_if(is.array, as.vector)

norm <- t(as.matrix(data[11:ncol(data)]))
colnames(norm) <- data$index

data_ann <- data[1:10]
rownames(data_ann) <- data$index

msce <- SingleCellExperiment(assays = list(logcounts = norm), 
                                          colData = data_ann)
rowData(msce)$feature_symbol <- rownames(msce)
colData(msce)$cell_type1 <- as.factor(data$louvain)
msce

geneIndex <- buildCellIndex(msce)
res <- findCell(geneIndex, genelist = gene_list)
res$common_exprs_cells


