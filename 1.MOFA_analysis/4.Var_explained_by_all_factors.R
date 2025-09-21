################################################################
# Plotting total variance explained by factor (all 50 factors) #
################################################################
library(MOFA2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(reshape2)
library(ggsci)

# Retrieving MOFA
load("C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Raw_data/mofas/standart_50f.RDATA")

# MOFA function to plot variance explained by factor
heat_1 <- plot_variance_explained(mofa, max_r2=15)$data

# Summarize, clean, sort, and format data for plotting
total_var_50f <- heat_1 %>%
  group_by(factor) %>%
  summarise(total_value = sum(value)) %>%
  mutate(factor = gsub('Factor','F',factor)) %>%
  arrange(desc(total_value)) %>%
  mutate(factor = factor(factor, level = unique(factor)))

factor_levels <- levels(total_var_50f$factor)

# Setting colors
color_gradient <- colorRampPalette(c("black", "blue"))(nrow(total_var_50f))

# Bar plot
xx = ggplot(total_var_50f,aes(x = factor, y=total_value, fill = factor))+
  geom_bar(stat = 'identity')+
  labs(x = 'Factors', y = '% Variance Explained')+
  scale_fill_manual(values = setNames(color_gradient, factor_levels))+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.text.x = element_text(color = ifelse(levels(total_var_50f$factor) %in% c('F1','F2','F3','F4','F8','F12','F14','F26','F42'), 'red', 'black'),
                                   size = 65, angle = 90),
        axis.title.x = element_text(size = 76),
        axis.title.y = element_text(size = 76),
        axis.text.y = element_text(size = 65))


# pdf(file = "C:/Users/User/Documents/Drive_PC/Lucas/GitHub/ROSMAP-MOFA/MOFA_real/Figuras_MOFA/pre_inkscape/geral/var_exp_factors_50.pdf",
#     width = 48, height = 27.5)
# xx
# dev.off()
