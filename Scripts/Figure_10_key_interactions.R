# ============================================================
# Figure 10: Key ligand-receptor interactions
# Final paper version summarizing the four highlighted CellChat interactions.
# ============================================================

library(CellChat)
library(dplyr)
library(ggplot2)

project_dir <- "C:/Ovarian_CAR_T_Project"
cellchat_path <- file.path(project_dir, "Figure_11_CellChat_revised_20261011", "checkpoint_sensitivity_truncatedMean_nboot100.rds")
figure_dir <- file.path(project_dir, "figures")
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

if (!file.exists(cellchat_path)) stop("CellChat object not found: ", cellchat_path)
cellchat <- readRDS(cellchat_path)
all_interactions <- subsetCommunication(cellchat, slot.name = "net", thresh = 1)
key_interactions <- all_interactions |>
  filter((ligand == "VEGFA" & receptor == "FLT1") | (ligand == "CXCL12" & receptor == "CXCR4") | (ligand == "THBS1" & receptor == "CD47") | (ligand == "THBS2" & receptor == "CD47")) |>
  mutate(interaction_label = paste0(ligand, "–", receptor), source_target = paste(source, "→", target)) |>
  arrange(desc(prob), interaction_label)
if (nrow(key_interactions) == 0) stop("No key interactions were found in the CellChat object.")
p <- ggplot(key_interactions, aes(x = reorder(interaction_label, prob), y = prob)) +
  geom_col(fill = "#4C78A8", width = 0.7) +
  geom_text(aes(label = source_target), hjust = -0.05, size = 3.3) +
  coord_flip() +
  labs(title = "Four key ligand–receptor interactions", x = NULL, y = "Communication probability") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(figure_dir, "Figure10_key_interactions.png"), p, width = 8.5, height = 5.5, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(figure_dir, "Figure10_key_interactions.pdf"), p, width = 8.5, height = 5.5, bg = "white", limitsize = FALSE)
write.csv(key_interactions, file.path(figure_dir, "Figure10_Four_Key_Interactions_truncatedMean.csv"), row.names = FALSE)
