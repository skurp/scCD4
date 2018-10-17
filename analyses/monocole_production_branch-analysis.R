# Will Connell
# 2018.09.24

library(monocle)
library(dplyr)


# Set environment ---------------------------------------------------------
# diretories
start.date <- Sys.Date()
analysis.date <- c("2018-09-25")
out.dir <- sprintf("../out/%s_monocle/%s_temporal-analysis", analysis.date, start.date)
dir.create(out.dir)
core_set <- round(detectCores()*0.5)

# Differential Expression Analysis 
# Cluster genes by pseudotemporal expression pattern
in.path <- sprintf("../out/%s_monocle/analyzed_obj.RDS", analysis.date)
HSMM <- readRDS(in.path)

diff_test_res <- differentialGeneTest(HSMM,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 1e-4))

png(sprintf("%s/pseudotemporal_expression.png", out.dir), width = 8, height = 18, units = 'in', res = 200)
plot_pseudotime_heatmap(HSMM,
                        num_clusters = length(unique(pData(HSMM)$Cluster)),
                        cores = core_set,
                        show_rownames = T)
dev.off()


# Branch point trajectory analysis
branch_traj_plots <- function(nbranches, HSMM) {
  for (branch in nbranches) {
    print(branch)
    BEAM_res <- BEAM(HSMM, branch_point = branch, cores = core_set)
    BEAM_res <- BEAM_res[order(BEAM_res$qval),]
    BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
    
    png(sprintf("%s/%i-branch.png", out.dir, branch), width = 8, height = 18, units = 'in', res = 200)
    plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res, qval < 1e-4)),],
                                branch_point = branch,
                                num_clusters = length(unique(pData(HSMM)$Cluster)),
                                cores = core_set,
                                use_gene_short_name = T,
                                show_rownames = T,
                                return_heatmap = F)
    dev.off()
  }
}

nbranches <- as.integer(sort(unique(pData(HSMM)$State)))
branch_traj_plots(nbranches, HSMM)


print('Done!')


# FUTURE GRANULAR ANALYSIS OF GENES
# test_genes <- row.names(subset(fData(HSMM),
#                                gene_short_name %in% c("CA5B", "KBTBD3", "ZNF605")))
# plot_genes_branched_pseudotime(HSMM[test_genes,],
#                                branch_point = 1,
#                                color_by = "louvain",
#                                ncol = 1)
