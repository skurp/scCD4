# Will Connell
# 2018.11.07

# Monocle 3 Setup

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("monocle")
devtools::install_github("cole-trapnell-lab/DDRTree", ref="simple-ppt-like")
devtools::install_github("cole-trapnell-lab/L1-graph")
install.packages("reticulate")
library(reticulate)
conda_create("monocle3")
use_condaenv("monocle3")
conda_install("monocle3", "umap-learn")
conda_install("monocle3", "louvain")
devtools::install_github("cole-trapnell-lab/monocle-release", ref="monocle3_alpha",
                         force = TRUE)
# test
library(monocle)



