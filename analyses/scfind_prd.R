# Will Connell
# 2018.09.28

# DE Genes across Clusters Interrogation

library(readr)
library(rhdf5)
library(scfind)
library(SingleCellExperiment)
library(tictoc)
library(dplyr)

# Prepare output paths ----------------------------------------------------
start.date <- Sys.Date()
out.dir <- sprintf("../out/%s_scfind-mat2", start.date)
dir.create(out.dir)

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
obs <- h5read(file_loc, "/obs") # observations, or rows, including metadata
vars <- h5read(file_loc, "/var") # variables, or columns

# subset cells by cluster indices
clust_indices <- which(obs$louvain %in% c(0,1))
# subset a sample to effectively
clust_indices <- sample(clust_indices, 20000)
h5closeAll()

print("Dimensions of FULL data:")
data <- h5read(file_loc, "/X") #normalized matrix
h5closeAll()

data <- data %>%
  t() %>% 
  as.data.frame()
colnames(data) <- c(vars$index)
data <- cbind(obs, data) %>% 
  mutate_if(is.array, as.vector)

print(str(data))

# ------- take subset of data to test function
data <- h5read(file_loc, "/X", index = list(NULL, clust_indices))
obs_indices <- obs[clust_indices,]

data <- data %>%
  t() %>% 
  as.data.frame()
colnames(data) <- c(vars$index)
data <- cbind(obs_indices, data) %>% 
  mutate_if(is.array, as.vector)

print("Dimensions of SUBSET data:")
print(str(data))

# ------ read in DE gene lists ---------
gene_loc <- "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.igtb.louvain.de.seurate.meta.louvain.sigpos.txt"
gene_data <- read_tsv(gene_loc)

print("Unique DE gene no. in each cluster:")
gene_data %>%
  count(cluster)
print(str(gene_data))

toc()


# SCfind setup ------------------------------------------------------------

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

find_cell_DE.genes <- function(geneIndex, clusters) {
  for (i in clusters) {
    gene_list <- gene_data %>% 
      filter(cluster == i) %>% 
      select(gene)
    res <- findCell(geneIndex, genelist = gene_list)
    png(sprintf("%s/%i-cluster_scfind_DE-search", out.dir, i), width = 8, height = 8, units = 'in', res = 200)
    barplot(-log10(res$p_values), ylab = "-log10(pval)", main = sprintf("HCA Gene Search: Cluster %i", i), las = 2)
    dev.off()
  }
}

find_cell_DE.genes(geneIndex, 0:1)

