# ============================================================
# Figure 11: Outgoing signaling patterns
# Final paper version from the CellChat signaling-role heatmap.
# ============================================================

library(CellChat)

project_dir <- "C:/Ovarian_CAR_T_Project"
cellchat_path <- file.path(project_dir, "Figure_11_CellChat_revised_20261011", "checkpoint_sensitivity_truncatedMean_nboot100.rds")
figure_dir <- file.path(project_dir, "figures")
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

if (!file.exists(cellchat_path)) stop("CellChat object not found: ", cellchat_path)
cellchat <- readRDS(cellchat_path)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

png(file.path(figure_dir, "Figure11_outgoing_signaling_heatmap.png"), width = 8.5, height = 7, units = "in", res = 300, bg = "white")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()

pdf(file.path(figure_dir, "Figure11_outgoing_signaling_heatmap.pdf"), width = 8.5, height = 7, bg = "white")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()
