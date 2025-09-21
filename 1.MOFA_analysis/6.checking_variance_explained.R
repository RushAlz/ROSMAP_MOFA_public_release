library(ggpubr)
library(MOFAdata)
library(MOFA2)
library(tidyverse)
setwd("/pastel/Github_scripts/ROSMAP-MOFA")

load("/pastel/projects/MOFA/model/standart_50f.RDATA")
phenotypes = readRDS("/pastel/resources/data2share/basic_Apr2022_selected_list_Jul2024.rds") 
phenotypes$projid = rownames(phenotypes)
load("/pastel/resources/data2share/pheno_list_Jul2024.RData")

sample_metadata <- data.frame(
  projid = samples_names(mofa)[[1]]
) %>% left_join(phenotypes)
sample_metadata$sample = sample_metadata$projid
samples_metadata(mofa) <- sample_metadata

loadings_scaled_ALL = get_weights(mofa, as.data.frame = T, scale = T, abs = T)

p_a = plot_variance_explained(mofa)

variance_explained_mtx = get_variance_explained(mofa)$r2_per_factor$group1
# scale each row (SUM = 1)
variance_explained_mtx_scaled = variance_explained_mtx / rowSums(variance_explained_mtx)
rowSums(variance_explained_mtx_scaled)

# For each factor (row) get the number of views (columns) that explain more than 0.1 of the variance (scaled) 
data.frame("n views" = rowSums(variance_explained_mtx_scaled > 0.1)) %>% 
  rownames_to_column("factor") -> n_views_per_factor
n_views_per_factor[n_views_per_factor$n.views==1,]

p_a2 = p_a + 
  geom_text(data = subset(p_a[[1]],value < 6), 
            aes(x = view, y = factor, 
                label = scales::number(value,accuracy = 0.1))) +
  geom_text(data = subset(p_a[[1]],value >= 6), 
            aes(x = view, y = factor, 
                label = scales::number(value,accuracy = 1)), color = "white") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "", x = "View", y = "Factor")

# Plot geom_bar
p_b = ggplot(n_views_per_factor, 
             aes(y = factor(factor, levels = n_views_per_factor$factor), x = n.views)) + geom_bar(stat = "identity") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = n.views), hjust = -0.5) +
  theme_classic() +
  labs(subtitle = "", x = "Num. views (>10% scaled)", y = "Factor")

pdf("figures/plot_variance_explained.pdf", width = 10, height = 12)
ggarrange(p_a2, p_b, ncol = 2, widths = c(2, 1))
dev.off()
