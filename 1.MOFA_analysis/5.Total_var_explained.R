##############################################
# Plotting total variance explained by views #
##############################################
library(MOFA2)
library(ggplot2)
library(dplyr)
library(ggsci)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# Retrieving participants phenotypes data
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/meta.RDATA")
meta <- meta[(meta$projid %in% mofa@samples_metadata$sample) , ]
samples_metadata(mofa) <- meta

# MOFA function to plot total variance explained
barrgraph_data <- plot_variance_explained(mofa, plot_total = T)[[2]]$data

# Reordering views
barrgraph_data <- barrgraph_data[,-3]
barrgraph_data <- barrgraph_data[c(7,1,2,3,5,6,4),]


# Renaming views
barrgraph_data$view <- c('H3K9ac', 'RNA\n(AC)','RNA\n(DLPFC)', 'RNA\n(PCG)', 'Protein', 
                         'Metabol.', 'Cell\nTypes')

barrgraph_data$view <- factor(barrgraph_data$view, levels = barrgraph_data$view)

# Bar colors
fill_col = pal_lancet(alpha = 1)(7)

# Bar plot
barrgraph <- ggplot(barrgraph_data, mapping = aes(x = view, y = R2)) +
  geom_col(fill = fill_col, color = 'black', width = 0.9) +
  labs(x = "Views", y = "% Variance", title = "Total Variance Explained by View") +
  theme_minimal() +  
  theme(plot.margin = margin(c(0,0,0,0)),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'plain'),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.grid.major.y= element_line(linewidth = 0.6),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 7, color = 'black'),
        axis.text.y = element_text(size = 7, color = 'black')) +
  scale_y_continuous(breaks = seq(0, 40, by = 10), expand = expansion(0,0)) +
  scale_x_discrete(expand = expansion(0.01,0))

barrgraph
