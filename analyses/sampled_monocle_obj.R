# Will Connell
# 2018.10.19

# Monocle Analysis - REPREX


# Depends -----------------------------------------------------------------
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("monocle")

library(rhdf5)
library(monocle)
library(dplyr)
library(readr)
library(tictoc)

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
clust_indices <- sample(1:nrow(obs), 20000)
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
sample_sheet <- data[1:10]
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



saveRDS(HSMM, "../data/sampled_monocle_obj.RDS")

print('DONE')
