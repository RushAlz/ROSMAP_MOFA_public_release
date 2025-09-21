#########################################
# Plotting variance explained by factor #
#########################################
library(MOFA2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(reshape2)
library(ggsci)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Retrieving participants phenotypes data
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/meta.RDATA")
meta <- meta[(meta$projid %in% mofa@samples_metadata$sample) , ]
samples_metadata(mofa) <- meta

# MOFA function to plot variance explained by factor
heat_1 <- plot_variance_explained(mofa, max_r2=15)$data

# Selecting AD-related factors only
heat_1 <- heat_1[which(heat_1$factor %in% c("Factor1", "Factor2", "Factor3", "Factor4", "Factor8", "Factor12", "Factor14", "Factor26", "Factor42")),]
heat_1 <- heat_1[,-4]

# Reshaping data.fame
heat_1 <- dcast(heat_1, factor ~ view, value.var = "value")
rownames(heat_1) <- heat_1$factor
heat_1 <- heat_1[,-1]

# Reordeing and renaming columns
heat_1 <- heat_1[,c(7,1,2,3,5,6,4)]
colnames(heat_1) <- c('H3K9ac', 'RNA (AC)', 'RNA (DLPFC)', 'RNA (PCG)', 'Proteomics',
                      'Metabolomics', 'Cell Types')

heat_1 <- as.matrix(heat_1)

# Setting heatmap colors
col <- pal_material('deep-orange')(10)
col <- c('white', col)
col <- colorRamp2(c(0,1,2,3,4,5,6,7,8,9,10), col)

# Heatmap
heat_var_xp <- Heatmap(heat_1,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "plain")),
                       name = '% of variance',
                       col = col,
                       column_title = "Variance Explained by Factor",
                       column_title_gp = gpar(fontsize = 12, fontface = 'plain'), 
                       rect_gp = gpar(col = 'black', lwd = 1, lty = 1),
                       border_gp = gpar(col = 'black', lwd=1),
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       row_order = rownames(heat_1)[order(as.integer(sub("Factor", "", rownames(heat_1))), decreasing = TRUE)],
                       column_order = paste(colnames(heat_1)),
                       row_names_side = 'left',
                       row_names_gp = gpar(fontsize = 9, fontface = 'plain'),
                       column_names_rot = 45,
                       column_names_gp = gpar(fontsize = 9),
                       heatmap_width = unit(60,'mm'),
                       heatmap_height = unit(75,'mm'))
print(heat_var_xp)

