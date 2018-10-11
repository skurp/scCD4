# Will Connell
# 2018.09.24

library(monocle)
library(dplyr)

# Differential Expression Analysis 

# Cluster genes by pseudotemporal expression pattern
HSMM <- readRDS("../out/analyzed_obj.RDS")

diff_test_res <- differentialGeneTest(HSMM,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 1e-4))

plot_pseudotime_heatmap(HSMM,
                        num_clusters = 9,
                        cores = 1,
                        show_rownames = T)

# Branch point trajectory analysis
BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 9,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

test_genes <- row.names(subset(fData(HSMM),
                               gene_short_name %in% c("CA5B", "KBTBD3", "ZNF605")))
plot_genes_branched_pseudotime(HSMM[test_genes,],
                               branch_point = 1,
                               color_by = "louvain",
                               ncol = 1)
